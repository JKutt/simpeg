from .....electromagnetics.static.resistivity.simulation import BaseDCSimulation as Sim
from .....utils import Zero, count, mkvc
from .....data import Data
from ....utils import compute_chunk_sizes
import warnings
from .....data import SyntheticData
import dask
import dask.array as da
import os
import shutil
import numpy as np
import time

Sim.sensitivity_path = './sensitivity/'


def dask_fields(self, m=None, compute_J=False):
    if m is not None:
        self.model = m
    strt = time.time()
    f = self.fieldsPair(self)
    # print('[info] fields pair: ', time.time() - strt)
    strt = time.time()
    A = self.getA()
    # print('[info] getA done: ', time.time() - strt)
    strt = time.time()
    Ainv = self.solver(A, **self.solver_opts)
    # print('[info] setup solver done: ', time.time() - strt)
    strt = time.time()
    RHS = self.getRHS()
    # print('[info] rhs done: ', time.time() - strt)
    strt = time.time()
    f[:, self._solutionType] = Ainv * RHS
    print('[info] solve done: ', time.time() - strt)
    strt = time.time()
    if compute_J:
        self.compute_J(f, Ainv)
    print('[info] J formed: ', time.time() - strt)
    return f


Sim.fields = dask_fields


def dask_dpred(self, m=None, f=None, compute_J=True):
    """
    dpred(m, f=None)
    Create the projected data from a model.
    The fields, f, (if provided) will be used for the predicted data
    instead of recalculating the fields (which may be expensive!).

    .. math::

        d_\\text{pred} = P(f(m))

    Where P is a projection of the fields onto the data space.
    """
    if self.survey is None:
        raise AttributeError(
            "The survey has not yet been set and is required to compute "
            "data. Please set the survey for the simulation: "
            "simulation.survey = survey"
        )

    if f is None:
        if m is None:
            m = self.model
        strt = time.time()
        f = self.fields(m, compute_J=compute_J)

    data = Data(self.survey)
    # strt = time.time()
    for src in self.survey.source_list:
        for rx in src.receiver_list:
            data[src, rx] = rx.eval(src, self.mesh, f)
    return mkvc(data)


Sim.dpred = dask_dpred


def dask_getJ(self, m, f=None):
    """
        Generate Full sensitivity matrix
    """
    if self._Jmatrix is None:
        self.fields(m, compute_J=True)

    return self._Jmatrix


Sim.getJ = dask_getJ


def dask_compute_J(self, f, Ainv):
    if os.path.exists(self.sensitivity_path):
        shutil.rmtree(self.sensitivity_path, ignore_errors=True)

        # Wait for the system to clear out the directory
        while os.path.exists(self.sensitivity_path):
            pass

    m_size = self.model.size
    count = 0
    dask_arrays = []
    for source in self.survey.source_list:
        u_source = f[source, self._solutionType]
        for rx in source.receiver_list:
            # wrt f, need possibility wrt m
            strt = time.time()
            PTv = rx.getP(self.mesh, rx.projGLoc(f)).toarray().T

            df_duTFun = getattr(f, "_{0!s}Deriv".format(rx.projField), None)
            df_duT, df_dmT = df_duTFun(source, None, PTv, adjoint=True)

            # Find a block of receivers
            n_block_col = int(np.ceil(df_duT.size * 8 * 1e-9 / self.max_ram))

            n_col = int(np.ceil(df_duT.shape[1] / n_block_col))

            nrows = int(
                m_size / np.ceil(m_size * n_col * 8 * 1e-6 / self.max_chunk_size)
            )
            # ind = 0
            # print('[info] getP done: ', time.time() - strt)
            strt = time.time()
            ATinvdf_duT = Ainv * df_duT
            # print('[info] ATinvdf_duT done: ', time.time() - strt)
            # strt = time.time()
            dA_dmT = self.getADeriv(u_source, ATinvdf_duT, adjoint=True)

            dRHS_dmT = self.getRHSDeriv(source, ATinvdf_duT, adjoint=True)

            du_dmT = da.from_delayed(
                    dask.delayed(-dA_dmT), shape=(m_size, n_col), dtype=float
            )
            if not isinstance(dRHS_dmT, Zero):
                du_dmT += da.from_delayed(
                    dask.delayed(dRHS_dmT), shape=(m_size, rx.nD), dtype=float
                )
            if not isinstance(df_dmT, Zero):
                du_dmT += da.from_delayed(
                    df_dmT, shape=(m_size, rx.nD), dtype=float
                )
            dask_arrays.append(du_dmT.T)
            del ATinvdf_duT
            # print('[info] blockin rx done: ', time.time() - strt)
            # strt = time.time()
            # for col in range(n_block_col):
            #     ATinvdf_duT = da.asarray(
            #         Ainv * df_duT[:, ind : ind + n_col]
            #     ).rechunk((nrows, n_col))

            #     dA_dmT = self.getADeriv(u_source, ATinvdf_duT, adjoint=True)

            #     dRHS_dmT = self.getRHSDeriv(source, ATinvdf_duT, adjoint=True)

            #     if n_col > 1:
            #         du_dmT = da.from_delayed(
            #             dask.delayed(-dA_dmT), shape=(m_size, n_col), dtype=float
            #         )
            #     else:
            #         du_dmT = da.from_delayed(
            #             dask.delayed(-dA_dmT), shape=(m_size,), dtype=float
            #         )

            #     if not isinstance(dRHS_dmT, Zero):
            #         du_dmT += da.from_delayed(
            #             dask.delayed(dRHS_dmT), shape=(m_size, n_col), dtype=float
            #         )

            #     if not isinstance(df_dmT, Zero):
            #         du_dmT += da.from_delayed(
            #             df_dmT, shape=(m_size, n_col), dtype=float
            #         )

            #     # blockName = self.sensitivity_path + "J" + str(count) + ".zarr"
            #     # da.to_zarr((du_dmT.T).rechunk("auto"), blockName)
            #     dask_arrays.append(du_dmT.T)
            #     del ATinvdf_duT
            #     count += 1

            #     ind += n_col
            # print('[info] blockin rx done: ', time.time() - strt)
    # dask_arrays = []
    # for ii in range(count):
    #     blockName = self.sensitivity_path + "J" + str(ii) + ".zarr"
    #     J = da.from_zarr(blockName)
    #     # Stack all the source blocks in one big zarr
    #     dask_arrays.append(J)

    # rowChunk, colChunk = compute_chunk_sizes(
    #     self.survey.nD, m_size, self.max_chunk_size
    # )
    # self._Jmatrix = da.vstack(dask_arrays).rechunk((rowChunk, colChunk))
    self._Jmatrix = da.vstack(dask_arrays)
    Ainv.clean()

    return self._Jmatrix


Sim.compute_J = dask_compute_J


def dask_getJtJdiag(self, m, W=None):
    """
        Return the diagonal of JtJ
    """
    if self.gtgdiag is None:

        # Need to check if multiplying weights makes sense
        if W is None:
            self.gtgdiag = da.sum(self.getJ(m) ** 2, axis=0).compute()
        else:
            w = da.from_array(W.diagonal())[:, None]
            self.gtgdiag = da.sum((w * self.getJ(m)) ** 2, axis=0).compute()

    return self.gtgdiag


Sim.getJtJdiag = dask_getJtJdiag


def dask_Jvec(self, m, v):
    """
        Compute sensitivity matrix (J) and vector (v) product.
    """
    self.model = m
    if self._Jmatrix is None:
        self.fields(m, compute_J=True)

    return da.dot(self._Jmatrix, v)


Sim.Jvec = dask_Jvec


def dask_Jtvec(self, m, v):
    """
        Compute adjoint sensitivity matrix (J^T) and vector (v) product.
    """
    self.model = m
    if self._Jmatrix is None:
        self.fields(m, compute_J=True)

    return da.dot(v, self._Jmatrix)


Sim.Jtvec = dask_Jtvec
