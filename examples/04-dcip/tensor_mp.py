import discretize
import discretize.utils as meshutils
from discretize.tensor_mesh import TensorMesh
from discretize.tree_mesh import TreeMesh
# from SimPEG import dask
# import dask
# from memory_profiler import profile
from multiprocessing import Process
import tracemalloc

from SimPEG import (
    maps,
    utils,
    data_misfit,
    regularization,
    optimization,
    inverse_problem,
    directives,
    inversion,
    objective_function,
    data
)
from SimPEG.electromagnetics.static import resistivity as dc, utils as DCutils
from dask.distributed import Client, LocalCluster
from dask import config
from SimPEG.utils.drivers import create_tile_meshes, create_nested_mesh
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from pymatsolver.direct import Pardiso as Solver
from pyMKL import mkl_set_num_threads
# from SimPEG.utils.solver_utils import SolverCHOLMOD as Solver
import time
import DCIPtools as DCIP
import os

def create_tile_dc(source, obs, uncert, global_mesh, global_active, tile_id):
    local_survey = dc.Survey(source)
    electrodes = np.vstack((local_survey.locations_a,
                            local_survey.locations_b,
                            local_survey.locations_m,
                            local_survey.locations_n))
    local_survey.dobs = obs
    local_survey.std = uncert

    # Create tile map between global and local
#     print('[info] creating nested mesh ')
#     time_nest = time.time()
    # local_mesh = create_nested_mesh(
    #     electrodes, global_mesh,
    # )
#     time_nest = time.time() - time_nest
#     print('[info] ', time_nest)

#     print('[info] time map')
#     time_map = time.time()
    # local_map = maps.TileMap(global_mesh, global_active, local_mesh)
    actmap = maps.InjectActiveCells(
        global_mesh, indActive=global_active, valInactive=np.log(1e-8)
    )
#     time_map = time.time() - time_map
#     print('[info] ', time_map)
    expmap = maps.ExpMap(global_mesh)
    mapping = expmap * actmap
    # Create the local misfit
    max_chunk_size = 256
    simulation = dc.Simulation3DNodal(
        global_mesh, survey=local_survey, sigmaMap=mapping, storeJ=False,
        Solver=Solver,
#         chunk_format="row",
#         max_chunk_size=max_chunk_size,
#         workers=workers
    )
    simulation.sensitivity_path = './sensitivity/Tile' + str(tile_id) + '/'
#     print(simulation.getSourceTerm().shape)
    data_object = data.Data(
        local_survey,
        dobs=obs,
        standard_deviation=uncert,
    )
    data_object.dobs = obs
    data_object.standard_deviation = uncert
    local_misfit = data_misfit.L2DataMisfit(
        data=data_object, simulation=simulation
    )
    local_misfit.W = 1 / uncert

    return local_misfit

# @profile
def run(start=0, end=0, cpumask={0, 1}):
    # mkl_set_num_threads(2)
    # os.fork()
    # client = Client(processes=True)
    # Get the set of CPUs
    # pid = os.getpid()
    # print("Number of CPUs:", os.cpu_count())
    # print("CPU affinity mask is for process id: ", pid, os.sched_getaffinity(pid))
    
    # os.sched_setaffinity(pid, cpumask)

    # print("CPU affinity mask is modified for process id: ", pid, os.sched_getaffinity(pid))
    

    out_prefix = "DC_4em3_40pct"
    h = [25, 25, 25]
    refine_topo = [0, 0, 4]
    refine_obs = [20, 10, 8, 8]
    refine_pad_obs = [10, 10, 4, 4]
    padDist = np.ones((3, 2)) * 1000
    rho_0 = 400


    # data_obj = DCutils.readUBC_DC3Dobs("Teck_Final_Decimated2_3pcUncer+4e-3floor_DC.obs")
    ct = time.time()
    patch = DCIP.loadDias(r"Teck_Final_All4.DAT")
    patch.readings = patch.readings[start:end]
    # patch.readings = patch.readings[0::30]
    print("[info] num readings: ", len(patch.readings))
    survey_dc = patch.createDcSurvey("DC", no_remote=False)
    data_obj = data.Data(survey_dc, dobs=survey_dc.dobs, standard_deviation=survey_dc.std)

    topo = np.loadtxt(r"LiDAR_points2.csv", delimiter=",")

    electrodes = utils.uniqueRows(np.vstack((data_obj.survey.locations_a,
                                data_obj.survey.locations_m,
                                data_obj.survey.locations_n)))
    electrodes = electrodes[0]

    # Calculate Geometric Factor
    G = DCutils.geometric_factor(data_obj.survey)
    # rhoApp = np.abs(data_obj.dobs * (1.0 / (G + 1e-10)))
    Vp = np.abs(rho_0 * G / np.abs(data_obj.dobs))

    # Reset uncerts
    data_obj.std = (np.abs(data_obj.dobs) + Vp) * 0.05
    print(f"New floor: {np.percentile(np.abs(data_obj.dobs), 10)}")


    # Quick drape
    tree = cKDTree(topo[:, :2])
    rad, ind = tree.query(electrodes[:, :2])
    electrodes[:, 2] = topo[ind, 2]

    b_pole = np.vstack(data_obj.survey.locations_b)
    rad, ind = tree.query(b_pole[:, :2])
    b_pole[:, 2] = topo[ind, 2]

    print("Creating global mesh")
    # strt_mesh = time()
    global_mesh = discretize.TensorMesh.readUBC(r"3DMesh_Fine.msh")

    print('[info] start active cells')
    active_cells = meshutils.active_from_xyz(global_mesh, topo)
    mstart = np.log(np.ones(int(active_cells.sum())) * (1.0 / rho_0))
    print('[info] draping electrodes ', global_mesh.nC, active_cells.sum())
    data_obj.survey.drape_electrodes_on_topography(global_mesh, active_cells, option='top')
    expmap = maps.ExpMap(global_mesh)
    mapactive = maps.InjectActiveCells(mesh=global_mesh, indActive=active_cells, valInactive=np.log(1e-8))
    mapping = expmap * mapactive
    global_mesh.writeUBC(f"{out_prefix}.msh", models={f"{out_prefix}_active.con": active_cells})
    print(f"Total number of sources: {len(data_obj.survey.source_list)}")
    print(f"Total number of data: {data_obj.survey.nD}")
    print('[info] start active cells info: ', global_mesh.nC, active_cells.sum())

    # ==============================================================================
    # finally split the simulations
    # -----------------------------

    print("[INFO] Creating tiled simulations over sources: ", len(survey_dc.source_list))
    eps = np.percentile(np.abs(survey_dc.dobs), 10)
    survey_dc.std = np.abs(survey_dc.dobs * 0.05) + eps
    survey = survey_dc
    # start_ = time.time()
    # local_misfits = []

    # idx_start = 0
    # idx_end = 0
    # # do every 5 sources
    # cnt = 0
    # src_collect = []
    # for ii, source in enumerate(survey.source_list):
    #     source._q = None # need this for things to work
    #     if cnt == 14 or ii == len(survey.source_list)-1:
    #         src_collect.append(source)        
    #         idx_end = idx_end + source.receiver_list[0].nD
    #         # dobs = survey_dc.dobs[idx_start:idx_end]
    # #         print(dobs.shape, len(src_collect))
    #         delayed_misfit = create_tile_dc(
    #                     src_collect,  survey.dobs[idx_start:idx_end],
    #                     survey.std[idx_start:idx_end], global_mesh, active_cells, ii)
    #         local_misfits += [delayed_misfit]
    # #         src_collect = client.scatter(src_collect)
    # #         survey_dobs = client.scatter(survey.dobs[idx_start:idx_end])
    # #         survey_std = client.scatter(survey.std[idx_start:idx_end])
    # #         global_mesh_ = client.scatter(global_mesh)
    # #         active_cells_ = client.scatter(active_cells)
    # #         ii = client.scatter(ii)
    # #         local_misfits += [client.submit(create_tile_dc,
    # #                     src_collect,  survey.dobs[idx_start:idx_end],
    # #                     survey.std[idx_start:idx_end], global_mesh, active_cells, ii)]
    #         idx_start = idx_end
    #         cnt = 0
    #         src_collect = []
    #     else:
    # #         print(idx_start, idx_end)
    #         src_collect.append(source)
    #         idx_end = idx_end + source.receiver_list[0].nD
    #         cnt += 1
    # # Simulations tiled map ==================================================
    # # local_misfits = client.gather(local_misfits)    

    # print("[INFO] Completed tiling in: ", time.time() - start_)
    # # local_misfits = local_misfits[:16]
    # electrodes_g = electrodes
    # global_misfit = objective_function.ComboObjectiveFunction(
    #                 local_misfits
    #         )
    # print("[info] number of misfit functions: ", len(local_misfits))
    simulation = dc.Simulation3DNodal(
        global_mesh, survey=survey_dc, sigmaMap=mapping, storeJ=False,
        Solver=Solver,)
    data_dc = data.Data(survey_dc, dobs=survey_dc.dobs, standard_deviation=survey_dc.std)
    global_misfit = data_misfit.L2DataMisfit(data=data_dc, simulation=simulation)
    # m0_dc = np.ones(active_cells.sum()) * np.log(1. / np.median(patch.getApparentResistivity()))
    # m0_dc = np.log(global_mesh.vol[active_cells])
    # Plot the model on different meshes
    ind = 6
    # fig = plt.figure(figsize=(14, 10))

    # local_mesh = global_misfit.objfcts[0].simulation.mesh
    # local_map = global_misfit.objfcts[0].simulation.sigmaMap
    # sub_survey = global_misfit.objfcts[0].simulation.survey
    # print(local_mesh.nC, global_mesh.nC)

    #====================================================================================
    # new implementation using the Combo Objective function and the dmis = dmis1 + dmis2
    # -----------------------------------------------------------------------------------
    #

    # make intital model
    m0_dc = np.ones(active_cells.sum()) * np.log(np.ones(int(active_cells.sum())) * (1.0 / rho_0))
    surface_weight = False
    use_preconditioner = False
    coolingFactor = 5
    coolingRate = 1
    beta0_ratio = 1e1

    # Map for a regularization
    regmap = maps.IdentityMap(nP=int(active_cells.sum()))
    # reg = regularization.Tikhonov(mesh, indActive=global_actinds, mapping=regmap)
    reg = regularization.Sparse(global_mesh, indActive=active_cells, mapping=regmap)

    # cell weights ==================================================
    if surface_weight:
        w_fac = np.ones(cell_depth_surface_w) * cell_weight
        surface_weights = get_surface_weights(global_mesh, active_cells, w_fac,
                                              octree=True)
        # Related to inversion
        reg.cell_weights = global_mesh.vol[active_cells] * surface_weights[active_cells]
        print('[INFO] Surface weights used')

    # print('[INFO] Getting things started on inversion...')
    # set alpha length scales
    reg.alpha_s = 1 # alpha_s
    reg.alpha_x = 1
    reg.alpha_y = 1
    reg.alpha_z = 1

    opt = optimization.ProjectedGNCG(maxIter=10, upper=np.inf, lower=-np.inf)
    invProb = inverse_problem.BaseInvProblem(global_misfit, reg, opt)

    # strt = time.time()
    # J = invProb.formJ(m0_dc)
    # # F = invProb.getFields(m0_dc)
    # print('time for fields: ', time.time() - strt)
    # # strt_dp = time.time()
    # # dpredd = invProb.get_dpred(m0_dc, f=None)
    # # # dmis_deriv = invProb.dmisfit.deriv(m0_dc, f=None)
    # # # print(invProb.dmisfit.objfcts[0].deriv(m0_dc, f=None))
    # # print('time for dpred: ', time.time() - strt_dp)

    # strt_d = time.time()
    # dmis_deriv = invProb.dmisfit.deriv(m0_dc, f=None)
    # print('time for deriv: ', time.time() - strt_d)

    beta = directives.BetaSchedule(
        coolingFactor=coolingFactor, coolingRate=coolingRate
    )
    betaest = directives.BetaEstimate_ByEig(beta0_ratio=beta0_ratio)
    target = directives.TargetMisfit()
    target.target = survey.nD / 2.
    saveIter = directives.SaveModelEveryIteration()
    saveIterVar = directives.SaveOutputEveryIteration()
    # Need to have basice saving function
    if use_preconditioner:
        update_Jacobi = directives.UpdatePreconditioner()
        updateSensW = directives.UpdateSensitivityWeights()
        # updateWj = Directives.Update_Wj()
        # directiveList = [
        #     beta, betaest, target, updateSensW, saveIter, update_Jacobi
        # ]
        directiveList = [
            updateSensW, beta, betaest, target, saveIter, update_Jacobi, saveIterVar 
        ]
    else:
        directiveList = [
            beta, betaest, target, saveIter
        ]
    inv = inversion.BaseInversion(
        invProb, directiveList=directiveList)
    opt.LSshorten = 0.5
    opt.remember('xc')
    print("[info] calc dpreds/fields")
    
    tracemalloc.start()
    return_g=True
    return_H=True
    outputs = invProb.getFields(m0_dc, store=(return_g is False and return_H is False))
    # invProb.phi_d = np.nan
    # invProb.phi_m = np.nan
    # outputs = invProb.evalFunction(m0_dc)
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()
    
    # print(outputs)
    

    # Run Inversion ================================================================
    # print("[INFO] starting inversion")
    # minv = inv.run(m0_dc)
    # print('time: ', time.time() - time_dpred)
    # rho_est = mapactive * minv
    # np.save('model_out.npy', rho_est)

if __name__ == '__main__':
    survey_type = 'dipole-dipole'
    # run(start=0, end=29)

    try:
        os.remove("flagfile.txt")
        os.remove("factorization_done.txt")
    except OSError:
        pass

    time_dpred = time.time()

    p1 = Process(target=run, kwargs={'start':0, 'end':300, 'cpumask':{0, 1, 2, 3, 4, 5}})
    # p2 = Process(target=run, kwargs={'start':100, 'end':200, 'cpumask':{6, 7, 8, 9, 10, 11}})
    # p3 = Process(target=run, kwargs={'start':200, 'end':300, 'cpumask':{8, 9, 10, 11}})

    p1.start()
    # p2.start()
    # p3.start()
    p1.join()
    # p2.join()
    # p3.join()

    print('time: ', time.time() - time_dpred)
