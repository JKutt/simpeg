#################################################
# Imports
import numpy as np
# from scipy import sparse
import matplotlib.pyplot as plt
from SimPEG.electromagnetics.static import resistivity as DC
from SimPEG.electromagnetics.static import induced_polarization as IP
from numpy import genfromtxt
# from openpyxl import load_workbook
from datetime import datetime
##################################################
# DCIPtools Version 2.2  Aug 2020

def loadDDN(file_input=None, rec_column='A', current_column='E', stn_column='C'):
    if file_input is not None:
        wb = load_workbook(filename=file_input)
        sheet = wb['InjectionsLog']
        print(sheet.max_row)
        # initiate list of active rows
        mem_i = []
        # open file for printing
        out_file = open(out_path, "w+")
        # lets find active rows for data extraction
        for idx in range(sheet.max_row):
            cell = rec_column + str(idx + 1)
            cell_in = current_column + str(idx + 1)
            cell_stn = stn_column + str(idx + 1)
            if isinstance(sheet[cell].value, int):
                if isinstance(sheet[cell_in].value, int):
                    mem_i.append([sheet[cell].value,
                                  sheet[cell_stn].value,
                                  sheet[cell_in].value])
                    # out_file.write('%i %s %i\n' % (sheet[cell].value,
                    #                                str(sheet[cell_stn].value),
                    #                                sheet[cell_in].value))
        # out_file.close()
        return mem_i
        print("done-zo!")
    else:
        print('[ERROR] Please supply an input DDN!')


def loadCsvDem(filepath):
    my_data = genfromtxt(filepath, delimiter=',', skip_header=1)
    utms = np.vstack((my_data[:, 2], my_data[:, 3]))
    elevations = my_data[:, 4]
    grid = np.vstack((my_data[:, 0], my_data[:, 1]))
    return utms, grid, elevations


def plotCsvDemLocationDifferences(local=None, utm=None):
    if local is None or utm is None:
        print("[INFO] Please provide Locations of DEM!!!!")
    else:
        delta_gridx = local[0, 1:] - local[0, :-1]
        delta_gridy = local[1, 1:] - local[1, :-1]
        delta_utmx = utm[0, 1:] - utm[0, :-1]
        delta_utmy = utm[1, 1:] - utm[1, :-1]

        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        ax1.plot(delta_gridx, '.')
        ax2.plot(delta_gridy, '.')
        ax3.plot(delta_utmx, '.')
        ax4.plot(delta_utmy, '.')

        ax1.set_title('Grid X difference')
        ax2.set_title('Grid Y difference')
        ax3.set_title('Easting difference')
        ax4.set_title('Northing difference')
        plt.show()


def removeRemoteFromGPS(gps_data, local, level):
    gps_edit = []
    local_edit = []
    level_edit = []
    for i in range(len(gps_data[0, :])):
        if local[0, i] == -9999 or local[0, i] == 9999:
            print("found remote")
        else:
            gps_edit.append(gps_data[:, i])
            local_edit.append(local[:, i])
            level_edit.append(level[i])
    return np.asarray(gps_edit).T, np.asarray(local_edit).T, np.asarray(level_edit).T


def plotLocalAndUtmGrids(local=None, utm=None):
    if local is None or utm is None:
        print("[INFO] Please provide locations info!!!")
    else:
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(local[0, :], local[1, :], '.')
        ax2.plot(utm[0, :], utm[1, :], '.')
        ax1.set_title('Local Grid Coords.')
        ax1.set_xlabel('X (m)')
        ax1.set_ylabel('Y (m)')
        ax2.set_xlabel('Easting (m)')
        ax2.set_ylabel('Northing (m)')
        ax2.set_title('UTM Coords.')
        ax2.legend(['local', 'utm'])
        plt.show()


def calcHomoHalfSpaceVoltage(background, tx1, tx2, rx1, rx2):
    r1 = ((rx1[0] - tx1[0])**2 +
          (rx1[1] - tx1[1])**2 +
          (rx1[2] - tx1[2])**2)**0.5
    r2 = ((rx2[0] - tx1[0])**2 +
          (rx2[1] - tx1[1])**2 +
          (rx2[2] - tx1[2])**2)**0.5
    r3 = ((rx1[0] - tx2[0])**2 +
          (rx1[1] - tx2[1])**2 +
          (rx1[2] - tx2[2])**2)**0.5
    r4 = ((rx2[0] - tx2[0])**2 +
          (rx2[1] - tx2[1])**2 +
          (rx2[2] - tx2[2])**2)**0.5
    In = 1                                           # normalized current
    v1 = (In * background) / (2 * np.pi * r1)        # voltage at point 1
    v2 = -1 * (In * background) / (2 * np.pi * r2)   # voltage at point 2
    v = v1 + v2                                      # theoretical voltage
    return v


def loadDias(fileName):
    """
    Function for loading a dias file and returns a
    "Patch" class complete with sources and recievers

    Input:
    fileName = complete path to data file

    """
    lines = 0
    text_file = open(fileName, "r")

    # determin how many lines in the file
    while text_file.readline():
            lines += 1
    text_file.close()

    # initiate a patch
    patch = Jpatch()
    # read header information
    text_file = open(fileName, "r")
    # initiate reading control variables
    currRdg = 0
    for i, line in enumerate(text_file):
        if i == 4:
            Varinfo = line.split()
            header5 = line
            # print(Varinfo)
        elif i == 0:
                header1 = line
        elif i == 1:
                header2 = line
                id_info = line.split()
                patch.assignPatchID(id_info[1])
        elif i == 2:
                header3 = line
        elif i == 3:
                header4 = line
        elif i > 4:
            try:
                datatxt = line.split()
                # do some Jdatamanagment stuff
                # print(i, line)
                varFields = Jreadtxtline(Varinfo, datatxt)
                # verify if line is a new reading
                if varFields.RDG == currRdg:
                    # add the dipoles
                    # Idp = Jdata.JinDipole(varFields)
                    Vpdp = JvoltDipole(varFields)
                    Rdg.addVoltageDipole(Vpdp)
                else:
                    # create a reading
                    Rdg = Jreading(varFields.RDG)
                    Idp = JinDipole(varFields)
                    Vpdp = JvoltDipole(varFields)
                    Rdg.addVoltageDipole(Vpdp)
                    Rdg.addInDipole(Idp)
                    # add reading to the patch
                    patch.addreading(Rdg)
                    currRdg = varFields.RDG
            except:
                    pass

    text_file.close()
    headers = [header1, header2, header3, header4, header5]
    patch.assignHeaderInfo(headers)
    return patch


# ===================================================
# Dias Data specific class
class JinDipole:
    """
    Class containing Source information

    Initiate with a structure containing location
    and Current value + error

    """

    def __init__(self, InDpInfo):
        self.Tx1File = InDpInfo.TxFile
        self.Tx1East = float(InDpInfo.Tx1East)
        self.Tx1North = float(InDpInfo.Tx1North)
        self.Tx1Elev = float(InDpInfo.Tx1Elev)
        self.Tx2East = float(InDpInfo.Tx2East)
        self.Tx2North = float(InDpInfo.Tx2North)
        self.Tx2Elev = float(InDpInfo.Tx2Elev)
        self.Tx1x = float(InDpInfo.Tx1x)
        self.Tx1y = float(InDpInfo.Tx1y)
        self.Tx2x = float(InDpInfo.Tx2x)
        self.Tx2y = float(InDpInfo.Tx2y)
        self.In = float(InDpInfo.In)
        self.In_err = float(InDpInfo.In_err)
        # print(InDpInfo.In_err)

    def getTxStack(self, stack_dir):
        rec_num = self.reading
        bin_file = stack_dir + "P1_R"+ rec_num + "_TX.raw"

        f = open(bin_file, "rb")          # open the file
        lines = f.read()                  # read data as string (can modify this to get header info
        return_char = 0.                  # beginning of binary data idex
        for idx in range(100):
            if lines[idx] == "\r":
                return_char = idx
        bin_start = return_char + 2       # increment 2 spaces to start of binary data
        data = []                         # initiate array for data
        with open(bin_file, 'rb') as f:   # open file for data extraction
            f.seek(bin_start,
                   os.SEEK_SET)          # seek to beginning of binary data
            while True:
                b = f.read(8)             # read 8 bytes at a time
                if not b:                 # break out if end of file
                    # eof
                    break
                data.append(
                    struct.unpack('d',
                                  b))  # store data

        return data


class JvoltDipole:
    """
    object containing voltage information

    """

    def __init__(self, VoltDpinfo):
        self.dipole = VoltDpinfo.DIPOLE
        self.reading = VoltDpinfo.RDG
        self.status = VoltDpinfo.Status
        self.Rx1x = float(VoltDpinfo.Rx1x)
        self.Rx1y = float(VoltDpinfo.Rx1y)
        self.Rx2x = float(VoltDpinfo.Rx2x)
        self.Rx2y = float(VoltDpinfo.Rx2y)
        self.Rx1File = VoltDpinfo.Rx1File
        self.Rx1East = float(VoltDpinfo.Rx1East)
        self.Rx1North = float(VoltDpinfo.Rx1North)
        self.Rx1Elev = float(VoltDpinfo.Rx1Elev)
        self.Rx2File = VoltDpinfo.Rx2File
        self.Rx2East = float(VoltDpinfo.Rx2East)
        self.Rx2North = float(VoltDpinfo.Rx2North)
        self.Rx2Elev = float(VoltDpinfo.Rx2Elev)
        self.coupling = float(VoltDpinfo.Coupling)
        self.K = float(VoltDpinfo.k)
        try:
            self.xcorr = float(VoltDpinfo.Xcorr)
        except:
            self.xcorr = -99.0
        try:
            self.eta = float(VoltDpinfo.M)
            self.c = float(VoltDpinfo.C)
            self.tau = float(VoltDpinfo.Tau)
        except:
            self.eta = -99.0
            self.c = -99.0
            self.tau = -99.0
        self.Sp = float(VoltDpinfo.Sp)
        try:
            self.Vp = float(VoltDpinfo.Vp)
        except:
            self.Vp = -999.9

        self.Vp_err = float(VoltDpinfo.Vp_err)
        self.Rho = float(VoltDpinfo.Rho)
        self.flagRho = VoltDpinfo.Rho_QC
        self.Stack = float(VoltDpinfo.Stack)
        try:
            self.Mx = float(VoltDpinfo.Mx)
        except:
            self.Mx = -99.9
        self.Mx_err = float(VoltDpinfo.Mx_err)
        self.flagMx = VoltDpinfo.Mx_QC
        self.flagBad = VoltDpinfo.Status
        self.TimeBase = VoltDpinfo.TimeBase
        self.Vs = np.asarray(VoltDpinfo.Vs)
        self.In = float(VoltDpinfo.In)
        self.In_err = float(VoltDpinfo.In_err)
        try:
            self.residual_dc = float(VoltDpinfo.rdc)
        except:
            self.residual_dc = -99.9
        try:
            self.residual_ip = float(VoltDpinfo.rip)
        except:
            self.residual_ip = -99.9

    def getDipoleStack(self, stack_dir):
        node1_id = self.Rx1File[:2]
        node2_id = self.Rx2File[:2]
        rec_num = self.reading
        bin_file = stack_dir + "P1_R"+ rec_num + "_" + node1_id + "_" + node2_id + ".stk"

        f = open(bin_file, "rb")          # open the file
        lines = f.read()                  # read data as string (can modify this to get header info
        return_char = 0.                  # beginning of binary data idex
        for idx in range(100):
            if lines[idx] == "\r":
                return_char = idx
        bin_start = return_char + 2       # increment 2 spaces to start of binary data
        data = []                         # initiate array for data
        with open(bin_file, 'rb') as f:   # open file for data extraction
            f.seek(bin_start,
                   os.SEEK_SET)           # seek to beginning of binary data
            while True:
                b = f.read(8)             # read 8 bytes at a time
                if not b:                 # break out if end of file
                    # eof
                    break
                data.append(
                    struct.unpack('d',
                                  b))  # store data
                # if len(b) < 8:
                #     break
        return data

    def getXplotpoint(self, Idp, dipole_dipole=False):
        x = -9999.
        if dipole_dipole:
            mid_tx = (Idp.Tx1x + Idp.Tx2x) / 2.0
            if (self.Rx1x > mid_tx):
                x = mid_tx + ((self.Rx1x - mid_tx) / 2.0)
            elif (self.Rx1x < mid_tx):
                x = Idp.Tx1x - ((Idp.Tx1x - self.Rx1x) / 2.0)
        else:
            if (self.Rx1x > Idp.Tx1x):
                x = Idp.Tx1x + ((self.Rx1x - Idp.Tx1x) / 2.0)
            elif (self.Rx1x < Idp.Tx1x):
                x = Idp.Tx1x - ((Idp.Tx1x - self.Rx1x) / 2.0)
        return x

    def getYplotpoint(self, Idp, dipole_dipole=False):
        y = -9999.
        if dipole_dipole:
            mid_tx = (Idp.Tx1y + Idp.Tx2y) / 2.0
            if (self.Rx1y > Idp.Tx1y):
                y = mid_tx + ((self.Rx1y - mid_tx) / 2.0)
            elif (self.Rx1y < mid_tx):
                y = mid_tx - ((mid_tx - self.Rx1y) / 2.0)
        else:
            if (self.Rx1y > Idp.Tx1y):
                y = Idp.Tx1y + ((self.Rx1y - Idp.Tx1y) / 2.0)
            elif (self.Rx1y < Idp.Tx1y):
                y = Idp.Tx1y - ((Idp.Tx1y - self.Rx1y) / 2.0)
        return y

    def getXplotpointUTM(self, Idp, dipole_dipole=False):
        x = -9999.
        if dipole_dipole:
            mid_tx = (Idp.Tx1East + Idp.Tx2East) / 2.0
            if (self.Rx1East > mid_tx):
                x = mid_tx + ((self.Rx1East - mid_tx) / 2.0)
            elif (self.Rx1East < mid_tx):
                x = Idp.Tx1East - ((Idp.Tx1East - self.Rx1East) / 2.0)
        else:
            if (self.Rx1East > Idp.Tx1East):
                x = Idp.Tx1East + ((self.Rx1East - Idp.Tx1East) / 2.0)
            elif (self.Rx1East < Idp.Tx1East):
                x = Idp.Tx1East - ((Idp.Tx1East - self.Rx1East) / 2.0)
        return x

    def getYplotpointUTM(self, Idp, dipole_dipole=False):
        y = -9999.
        if dipole_dipole:
            mid_tx = (Idp.Tx1North + Idp.Tx2North) / 2.0
            if (self.Rx1North > Idp.Tx1North):
                y = mid_tx + ((self.Rx1y - mid_tx) / 2.0)
            elif (self.Rx1North < mid_tx):
                y = mid_tx - ((mid_tx - self.Rx1North) / 2.0)
        else:
            if (self.Rx1North > Idp.Tx1North):
                y = Idp.Tx1North + ((self.Rx1North - Idp.Tx1North) / 2.0)
            elif (self.Rx1North < Idp.Tx1North):
                y = Idp.Tx1North - ((Idp.Tx1North - self.Rx1North) / 2.0)
        return y

    def getZplotpoint(self, Idp, dipole_dipole=False):
        z = -9999.
        if dipole_dipole:
            # mid_tx = (Idp.Tx1x + Idp.Tx2x) / 2.0
            # z = -(abs(mid_tx - self.Rx1x)) / 2.0
            # mid_tx = (Idp.Tx1y + Idp.Tx2y) / 2.0
            # z = -(abs(mid_tx - self.Rx1y)) / 2.0
            r4 = ((self.Rx2East - Idp.Tx2East)**2 +
              (self.Rx2North - Idp.Tx2North)**2 +
              (self.Rx2Elev - Idp.Tx2Elev)**2)**0.5
            z = r4/ 3.
        else:
            r = np.sqrt((Idp.Tx1East - self.Rx2East)**2 +
                        (Idp.Tx1North - self.Rx2North)**2 +
                        (Idp.Tx1Elev - self.Rx2Elev)**2)
            z = -(r / 3.)
        return z

    def calcSensitivity(self, Idp, ne, sw, dx,
                        sw_search,
                        ne_search):
        xa = Idp.Tx2East
        ya = Idp.Tx2North
        za = Idp.Tx2Elev
        xb = Idp.Tx1East
        yb = Idp.Tx1North
        zb = Idp.Tx1Elev
        xm = self.Rx1East
        ym = self.Rx1North
        zm = self.Rx1Elev
        xn = self.Rx2East
        yn = self.Rx2North
        zn = self.Rx2Elev
        X = np.arange(sw[0], ne[0], dx)
        Y = np.arange(sw[1], ne[1], dx)
        Z = np.arange(sw[2], ne[2], dx)
        x, y, z = np.meshgrid(X, Y, Z)
        # compactify the calculation
        arg1 = (((x - xa) * (x - xm) + (y - ya) * (y - ym) +
                (z - za) * (z - zm)) *
                (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
                (((x - xm)**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
        arg2 = (((x - xa) * (x - xn) + (y - ya) * (y - yn) +
                (z - za) * (z - zn)) *
                (((x - xa)**2 + (y - ya)**2 + (z - za)**2)**-1.5) *
                (((x - xn)**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
        arg3 = (((x - xb) * (x - xm) + (y - yb) * (y - ym) +
                (z - zb) * (z - zm)) *
                (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
                (((x - xm)**2 + (y - ym)**2 + (z - zm)**2)**-1.5))
        arg4 = (((x - xb) * (x - xn) + (y - yb) * (y - yn) +
                (z - zb) * (z - zn)) *
                (((x - xb)**2 + (y - yb)**2 + (z - zb)**2)**-1.5) *
                (((x - xn)**2 + (y - yn)**2 + (z - zn)**2)**-1.5))
        # evaluate sensitivity equation
        J = 1 / (4 * np.pi**2) * (arg1 - arg2 - arg3 + arg4)
        # get the span of the most sensitive layer
        max_J = np.max(np.abs(J)) * 0.91
        levels = J.shape[2]
        slices_x = np.zeros((levels, 2))
        slices_y = np.zeros((levels, 2))
        for l in range(levels):
            z = J[:, :, l]
            a = np.abs(z) > (np.max(np.abs(z)) * 0.61)
            X_ = x[a]
            Y_ = y[a]
            try:
                slices_x[l, :] = np.asarray([np.max(X_), np.min(X_)])
                slices_y[l, :] = np.asarray([np.max(Y_), np.min(Y_)])
            except ValueError:
                slices_x[l, :] = slices_x[0, :]
                slices_y[l, :] = slices_y[0, :]

        sbx = [np.max(slices_x[:, 0]), np.min(slices_x[:, 1])]
        sby = [np.max(slices_y[:, 0]), np.min(slices_y[:, 1])]
        # determine if it overlaps with the area in question
        true_x = False
        true_y = False
        x_s = sw_search[0]
        x1_s = ne_search[0]
        y_s = sw_search[1]
        y1_s = ne_search[1]
        if (sbx[1] < x_s < sbx[0]):
            true_x = True
        if (sbx[1] < x1_s < sbx[0]):
            true_x = True
        if (sby[1] < y_s < sby[0]):
            true_y = True
        if (sby[1] < y1_s < sby[0]):
            true_y = True
        # print("Max sensitive x comp. {0} - min: {1}".format(sbx[0], sbx[1]))
        # print("Max search x comp. {0} - min: {1}".format(x1_s, x_s))
        # print("Max sensitive y comp. {0} - min: {1}".format(sby[0], sby[1]))
        # print("Max search y comp. {0} - min: {1}".format(y1_s, y_s))
        # print(true_x, true_y)
        if true_x and true_y:
            print("found sensitive dipole")
            self.flagRho = "Reject"
            self.flagMx = "Reject"

        return J

    def getDepthPoint2D(self, Idp, direction='x', dipole_dipole=False):
        z = np.nan
        if dipole_dipole:
            if direction is 'x':
                # mid_tx = (Idp.Tx1x + Idp.Tx2x) / 2.
                # mid_rx = (self.Rx1x + self.Rx2x) / 2

                # z = -(abs(mid_tx - mid_rx)) / 2.0
                z = (Idp.Tx1x + Idp.Tx2x + self.Rx1x + self.Rx2x) / 4.
            elif direction is 'y':
                # mid_tx = (Idp.Tx1y + Idp.Tx2y) / 2.
                # mid_rx = (self.Rx1y + self.Rx2y) / 2

                # z = -(abs(mid_tx - mid_rx)) / 2.0
                z = (Idp.Tx1y + Idp.Tx2y + self.Rx1y + self.Rx2y) / 4.
        else:
            if direction is 'x':
                z = -(abs(Idp.Tx1x - self.Rx1x)) / 2.0
            elif direction is 'y':
                z = -(abs(Idp.Tx1y - self.Rx1y)) / 2.0
        return z

    def getNlevel(self, Idp, a, direction='x'):
        n = np.nan
        if direction is 'x':
            if Idp.Tx1x > self.Rx1x:
                if Idp.Tx1x < Idp.Tx2x:
                    n = (abs(Idp.Tx1x - self.Rx1x)) / a
                else:
                    n = (abs(Idp.Tx2x - self.Rx1x)) / a
            else:
                if Idp.Tx1x > Idp.Tx2x:
                    n = (abs(Idp.Tx1x - self.Rx1x)) / a
                else:
                    n = (abs(Idp.Tx2x - self.Rx1x)) / a
        elif direction is 'y':
            if Idp.Tx1y > self.Rx1y:
                if Idp.Tx1y < Idp.Tx2y:
                    n = (abs(Idp.Tx1y - self.Rx1y)) / a
                else:
                    n = (abs(Idp.Tx2y - self.Rx1y)) / a
            else:
                if Idp.Tx1y > Idp.Tx2y:
                    n = (abs(Idp.Tx1y - self.Rx1y)) / a
                else:
                    n = (abs(Idp.Tx2y - self.Rx1y)) / a
        return n

    def getAseperation(self):
        r1 = ((self.Rx1East - self.Rx2East)**2 +
              (self.Rx1North - self.Rx2North)**2 +
              (self.Rx1Elev - self.Rx2Elev)**2)**0.5
        return r1

    def getTxRxSeperation(self, Idp):
        r1 = ((self.Rx1East - Idp.Tx1East)**2 +
              (self.Rx1North - Idp.Tx1North)**2 +
              (self.Rx1Elev - Idp.Tx1Elev)**2)**0.5
        return r1

    def getAseperationLocal(self, direction='x'):
        a = np.nan
        if direction is 'x':
            a = np.abs(self.Rx2x - self.Rx1x)
        elif direction is 'y':
            a = np.abs(self.Rx2y - self.Rx1y)
        return a

    def calcAngleFromNorth(self, Idp):
        # alpha = self.calcAngleFromTx(Idp)
        x_diff = self.Rx2East - self.Rx1East
        y_diff = self.Rx2North - self.Rx1North

        gamma = (np.arctan2(x_diff, y_diff)) * 180 / np.pi
        if gamma < 0:
            gamma = 360 + gamma
        # if (alpha < 120) or (alpha > 300):
        #     if (gamma > 120) or (gamma < 300):
        #         gamma = gamma - 180
        # else:
        #     if (gamma < 120) or (gamma > 300):
        #         gamma = gamma + 180

        return gamma

    def calcAngleFromTx(self, Idp):
        midx = (self.Rx1East + self.Rx2East) / 2.0
        midy = (self.Rx1North + self.Rx2North) / 2.0

        x_diff = midx - Idp.Tx1East
        y_diff = midy - Idp.Tx1North

        # x_diff = self.Rx2East - Idp.Tx1East
        # y_diff = self.Rx2North - Idp.Tx1North

        alpha = (np.arctan2(x_diff, y_diff)) * 180 / np.pi
        if alpha < 0:
            alpha = 360 + alpha
        return alpha

    def calcGeoFactor(self, Idp):
        r1 = ((self.Rx1East - Idp.Tx1East)**2 +
              (self.Rx1North - Idp.Tx1North)**2 +
              (self.Rx1Elev - Idp.Tx1Elev)**2)**0.5
        r2 = ((self.Rx2East - Idp.Tx1East)**2 +
              (self.Rx2North - Idp.Tx1North)**2 +
              (self.Rx2Elev - Idp.Tx1Elev)**2)**0.5
        r3 = ((self.Rx1East - Idp.Tx2East)**2 +
              (self.Rx1North - Idp.Tx2North)**2 +
              (self.Rx1Elev - Idp.Tx2Elev)**2)**0.5
        r4 = ((self.Rx2East - Idp.Tx2East)**2 +
              (self.Rx2North - Idp.Tx2North)**2 +
              (self.Rx2Elev - Idp.Tx2Elev)**2)**0.5
        try:
            gf = 1 / ((1 / r1 - 1 / r2) - (1 / r3 - 1 / r4))
        except ZeroDivisionError:
            gf = 0.0
        return 2 * np.pi * gf

    def calcRho(self, Idp):
        r1 = ((self.Rx1East - Idp.Tx1East)**2 +
              (self.Rx1North - Idp.Tx1North)**2 +
              (self.Rx1Elev - Idp.Tx1Elev)**2)**0.5
        r2 = ((self.Rx2East - Idp.Tx1East)**2 +
              (self.Rx2North - Idp.Tx1North)**2 +
              (self.Rx2Elev - Idp.Tx1Elev)**2)**0.5
        r3 = ((self.Rx1East - Idp.Tx2East)**2 +
              (self.Rx1North - Idp.Tx2North)**2 +
              (self.Rx1Elev - Idp.Tx2Elev)**2)**0.5
        r4 = ((self.Rx2East - Idp.Tx2East)**2 +
              (self.Rx2North - Idp.Tx2North)**2 +
              (self.Rx2Elev - Idp.Tx2Elev)**2)**0.5
        try:
            gf = 1 / ((1 / r1 - 1 / r2) - (1 / r3 - 1 / r4))
        except ZeroDivisionError:
            gf = 0.0
        Vp = self.Vp
        # print("Vp: {0}".format(self.Vp))
        rho = (Vp / self.In) * 2 * np.pi * gf
        # self.Rho = rho
        return rho

    def calcMx(self, window_widths, start, stop):
        widths = window_widths[start:stop]
        vs_ = self.Vs[start:stop]
        mx_calc = np.sum(vs_ * widths) / np.sum(widths)
        if self.Vp != 0:
            mx_calc = mx_calc / (self.Vp / 1000.0)
        else:
            mx_calc = -99.0
        return mx_calc


class Jreading:
    """
    Class to handle current and voltage dipole
    information for a given source

    """

    def __init__(self, mem):
        self.MemNumber = mem
        self.Vdp = []
        self.win_width = []
        self.win_start = []
        self.win_end = []
        self.node_db = []
        self.node_locs = []
        self.node_ids = []
        self.angles = []

    # method for adding voltage dipoles
    def addVoltageDipole(self, JVdp):
        self.Vdp.append(JVdp)

    # method for assigning Current dipole
    def addInDipole(self, JIdp):
        self.Idp = JIdp

    # method for creating a node database
    def createNodeDB(self):
        nodes = []
        locs = []
        node_num = []
        for dp in range(len(self.Vdp)):
            node1_id = self.Vdp[dp].Rx1File[:2]
            node2_id = self.Vdp[dp].Rx2File[:2]
            if not node1_id in nodes:
                nodes.append(node1_id)
                pos = np.asarray([self.Vdp[dp].Rx1East, self.Vdp[dp].Rx1North])
                locs.append(pos)
                node_num.append(self.Vdp[dp].Rx1File)
            if not node2_id in nodes:
                nodes.append(node2_id)
                pos = np.asarray([self.Vdp[dp].Rx2East, self.Vdp[dp].Rx2North])
                locs.append(pos)
                node_num.append(self.Vdp[dp].Rx2File)
        # assign node database
        self.node_db = nodes
        self.node_locs = locs
        self.node_ids = node_num

    # method for getting angles of possible dipole given a spacing
    def calcAngles(self, min_dipole, max_dipole):
        # dipoles = [] len(self.node_locs
        if len(self.node_locs) > 1:
            points = np.asarray(self.node_locs)
            for idx in range(len(self.node_locs)):
                # calculate difference to all points
                diff = self.node_locs[idx] - points
                # calculate distances
                dists = (np.sum(diff**2, axis=1))**0.5
                for index in range(dists.size):
                    if min_dipole < dists[index] < max_dipole:
                        # get components differences
                        x_diff = diff[index, 0]
                        y_diff = diff[index, 1]
                        # calc angle from north
                        theta = (np.arctan2(x_diff, y_diff)) * 180 / np.pi
                        if theta < 0:
                            theta = 360 + theta
                        self.angles.append(theta)
        else:
            print("node locations have not been initiated")


class Jpatch:
    """
    Class to hold source information for a given data patch

    """

    def __init__(self):
        self.readings = []

    def addreading(self, Jrdg):
        self.readings.append(Jrdg)

    def assignPatchID(self, id):
        self.ID = id

    def assignHeaderInfo(self, headerLines):
        """
        Method for processing the header lines
        of a Dias data file. (e.g IP times, project name, etc.)

        Input: an array of header lines from file

        """

        self.headers = headerLines         # assigns headers to patch class
        # process IP times from header
        timing_string = self.headers[2].split(' ')[1]
        timings = timing_string.split(';')[0]
        timing = timings.split(',')
        self.window_start = np.zeros(len(timing))
        self.window_end = np.zeros(len(timing))
        self.window_center = np.zeros(len(timing))
        self.window_width = np.zeros(len(timing))
        # loops each available IP window and calculate width and centre
        for i in range(len(timing)):
            wintime = timing[i].split(':')
            self.window_start[i] = float(wintime[0])
            self.window_end[i] = float(wintime[1])
            self.window_center[i] = (self.window_start[i] +
                                     self.window_end[i]) / 2.0
            self.window_width[i] = (self.window_end[i] -
                                    self.window_start[i])

    def assignWindowWidthsToDipoles(self):
        """
        assigns widths of windows of Vs decay to each dipole

        """
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            self.readings[k].win_width = self.window_width
            self.readings[k].win_start = self.window_start
            self.readings[k].win_end = self.window_end

    # def get3DanglesFromTx(self):
    #     theta = []
    #     num_rdg = len(self.readings)
    #     for rdg in range(num_rdg):
    #         In = self.readings[rdg].Idp
    #         num_dipole = len(self.readings[rdg].Vdp)
    #         for dp in range(num_dipole):
    #             # get a node database for a reading
    #             alpha = self.readings[rdg].Vdp[dp].calcAngleFromTx(In)
    #             if alpha <= 0:
    #                 alpha = alpha + 360
    #             theta.append(alpha)
    #             self.readings[rdg].angles.append(alpha)
    #     return theta

    def get3DanglesFromTx(self, line_azimuth):
        theta = []
        azimuth1 = line_azimuth
        azimuth2 = line_azimuth + 180
        num_rdg = len(self.readings)
        for rdg in range(num_rdg):
            In = self.readings[rdg].Idp
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                # get a node database for a reading
                alpha = self.readings[rdg].Vdp[dp].calcAngleFromTx(In)
                alpha1 = self.readings[rdg].Vdp[dp].calcAngleFromNorth(In)

                if alpha1 < 0:
                    alpha1 = alpha1 + 360
                if azimuth1 < alpha < 180:
                    if azimuth1 > alpha1 > 0:
                        alpha1 = alpha1 + 180
                    if 360 > alpha1 < azimuth2:
                        alpha1 = alpha1 - 180
                elif azimuth2 > alpha > 180:
                    if 360 > alpha1 > azimuth2:
                        alpha1 = alpha1 - 180
                    if azimuth1 > alpha1 > 0:
                        alpha1 = alpha1 + 180
                elif azimuth1 > alpha > 0:
                    if 180 > alpha1 > azimuth1:
                        alpha1 = alpha1 + 180
                    if azimuth2 > alpha1 >= 180:
                        alpha1 = alpha1 - 180
                elif 360 > alpha > azimuth2:
                    if azimuth2 > alpha1 > 180:
                        alpha1 = alpha1 - 180
                    if 180 > alpha1 > azimuth1:
                        alpha1 = alpha1 + 180
                if alpha1 < 0:
                    alpha1 = alpha1 + 360
                theta.append(alpha1)
                self.readings[rdg].angles.append(alpha1)
        return theta

    def rejectByVsError(self):
        num_rdg = len(self.readings)
        for rdg in range(num_rdg):
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                if self.readings[rdg].Vdp[dp].Mx_err == -99.9:
                        self.readings[rdg].Vdp[dp].flagMx = "Reject"

    def get3Dangles(self, bin_width=20, dipole_max=None, dipole_min=None):
        dipole_angles = []
        theta = []
        cnt = 0
        if dipole_max is not None:
            num_rdg = len(self.readings)
            for rdg in range(num_rdg):
                # get a node database for a reading
                self.readings[rdg].createNodeDB()
                self.readings[rdg].calcAngles(dipole_min, dipole_max)
                dipole_angles.append(self.readings[rdg].angles)
                for idx in range(len(dipole_angles[cnt])):
                    theta.append(dipole_angles[cnt][idx])
                cnt = cnt + 1
        else:
            print("[OUTPUT] - dipole sizes not set")
        return theta

    # method for making node database for entire project
    def getNodeDataBase(self):
        num_rdg = len(self.readings)
        nodes = []
        # locs = []
        for rdg in range(num_rdg):
            self.readings[rdg].createNodeDB()
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                node1_id = self.readings[rdg].Vdp[dp].Rx1File[:2]
                node2_id = self.readings[rdg].Vdp[dp].Rx2File[:2]
                if not node1_id in nodes:
                    nodes.append(node1_id)
                    pos = np.asarray([self.readings[rdg].Vdp[dp].Rx1East,
                                     self.readings[rdg].Vdp[dp].Rx1North])
                    # locs.append(pos)
                if not node2_id in nodes:
                    nodes.append(node2_id)
                    pos = np.asarray([self.readings[rdg].Vdp[dp].Rx2East,
                                     self.readings[rdg].Vdp[dp].Rx2North])
                    # locs.append(pos)
        return nodes

    def checkForMissingRecords(self, record_log=None):
        if record_log is not None:
            rec_out_file = record_log + ".rlog"
            text_file = open(record_log, "r")
            outfile = open(rec_out_file, "w+")
            # retrieve and assign the data
            rec_log = []
            current_rec = -1
            for line in text_file:
                spl = line.split(':')
                if int(spl[1]) != current_rec:
                    rec_log.append(int(spl[1]))
                    current_rec = int(spl[1])
            for rdg in range(len(self.readings)):
                found_match = False
                for rec in rec_log:
                    if int(self.readings[rdg].MemNumber) == rec:
                        found_match = True
                if not found_match:
                    outfile.write('rec %s found no match\n' % self.readings[rdg].MemNumber)
            print("[OUTPUT]  Log file written to same file as records log location")
        else:
            print(      "Please supply a records log file!!!!!")

    def checkContinuityInNodeFiles(self, path=None):
        if path is not None:
            nodes = self.getNodeDataBase()
            nodes.sort()
            # print(nodes)
            outfile = open(path, "w+")
            for idx in range(len(nodes)):
                nid = nodes[idx]
                files = []
                hold_num = -99
                hold_id = ""
                hold_mem = -99
                for rdg in range(len(self.readings)):
                    for ids in range(len(self.readings[rdg].node_ids)):
                        if nid == self.readings[rdg].node_ids[ids][:2]:
                            if hold_num == -99:
                                outfile.write('node %s mem start: %s\n' % (nid, self.readings[rdg].MemNumber))
                                hold_num = int(self.readings[rdg].node_ids[ids][2:8])
                                hold_mem = int(self.readings[rdg].MemNumber)
                                hold_id = self.readings[rdg].node_ids[ids]
                            else:
                                diff = int(self.readings[rdg].node_ids[ids][2:8]) - hold_num
                                if diff > 1:
                                    outfile.write('found gap: %s mem: %i & %s mem: %s\n' % (hold_id,
                                                   hold_mem, self.readings[rdg].node_ids[ids],
                                                   self.readings[rdg].MemNumber))
                                    hold_num = int(self.readings[rdg].node_ids[ids][2:8])
                                    hold_mem = int(self.readings[rdg].MemNumber)
                                    hold_id = self.readings[rdg].node_ids[ids]
            # write node database collected from db file
            for ids in nodes:
                outfile.write('%s\n' % ids)
            print("It be done!, file written")
        else:
            print("please supply path to database file")

    def plotElevations(self, utm=None, line_direction=None, dx=10):
        rx_locs = self.getDipoles(reject=None)
        tx_locs = self.getSources(reject=None)
        sw = np.array([np.min(rx_locs[:, 0]), np.min(tx_locs[:, 1])])
        ne = np.array([np.max(rx_locs[:, 0]), np.max(rx_locs[:, 1])])
        X = np.arange(sw[0], ne[0], dx)
        Y = np.arange(sw[1], ne[1], dx)
        x, y = np.meshgrid(X, Y)
        z = np.zeros_like(x)

        # need to find elevations by creating array
        levels = []
        locs = []
        for rdg in range(len(self.readings)):
                current = self.readings[rdg].Idp
                print("rec num: {0}".format(self.readings[rdg].MemNumber))
                levels.append(self.readings[rdg].Idp.Tx1Elev)
                locs.append([self.readings[rdg].Idp.Tx1East, self.readings[rdg].Idp.Tx1North])
                for dp in range(len(self.readings[rdg].Vdp)):
                    levels.append(self.readings[rdg].Vdp[dp].Rx1Elev)
                    levels.append(self.readings[rdg].Vdp[dp].Rx2Elev)
                    locs.append([self.readings[rdg].Vdp[dp].Rx1East, self.readings[rdg].Vdp[dp].Rx1North])
                    locs.append([self.readings[rdg].Vdp[dp].Rx2East, self.readings[rdg].Vdp[dp].Rx2North])

    def rejectDataFromSensitiveArea(self, dx, sw_s, ne_s, outFile=None):
        """
           This function takes in an area of concern and removes data
           sensitive to the defined area
        """
        if outFile is not None:
            rx_locs = self.getDipoles()
            tx_locs = self.getSources()
            sw = np.array([np.min(rx_locs[:, 0]), np.min(tx_locs[:, 1]),
                           np.min(rx_locs[:, 2]) - 200])
            ne = np.array([np.max(rx_locs[:, 0]), np.max(rx_locs[:, 1]),
                           np.max(rx_locs[:, 2])])
            X = np.arange(sw[0], ne[0], dx)
            Y = np.arange(sw[1], ne[1], dx)
            Z = np.arange(sw[2], ne[2], dx)
            J_ = np.zeros((Y.size, X.size, Z.size))
            for rdg in range(len(self.readings)):
                current = self.readings[rdg].Idp
                print("rec num: {0}".format(self.readings[rdg].MemNumber))
                for dp in range(len(self.readings[rdg].Vdp)):
                    J_ += self.readings[rdg].Vdp[dp].calcSensitivity(current, ne, sw, dx, sw_search=sw_s, ne_search=ne_s)
            
            start_time = self.window_start[0]
            end_time = self.window_end[0]
            start_inds = (window_start >= start_time)
            stop_inds = (window_end <= end_time)
            strt_time = window_start[start_inds]
            stop_time = window_end[stop_inds]
            # write the data to file
            self.writeColeColeSEDat(outname, strt_time[0],
                                     stop_time[stop_time.size - 1])
            print("New File written with sensitive area removed")
        else:
            print("     Please supply an output file path!")       
        
    def writeOldFartDat(self, outname, start_time, end_time, direction, a, line):
        """
        Writes data to old man 2D data for DiasQC
        """
        now = datetime.now()
        date = now.strftime("%d/%m/%Y")

        projectName  = outname.split('/')[-1]

        broken = outname.split('.')
        out2 = broken[0] + "-flag.DAT"
        out_file = open(outname, "w+")
        format_jup = '%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10i %10i %10i %10i %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.2f %10.2f %12.3f %10.3f %10.4f %10.2f %10.2f %8i %10.2f'
        format_flag = '%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10i %10i %10i %10i %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.2f %10.2f %12.3f %10.3f %10.4f %10.2f %10.2f %8i %10s'
        format_dc = '%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10i %10i %10i %10i %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.2f %10.2f %12.3f %10.3f %10.4f %10s %10s %8i %10s'
        # write the headers
        header_columns = "       C1X        C1Y        C2X        C2Y        P1X        P1Y        P2X        P2Y    Tx1Line    Tx2Line    Rx1Line    Rx2Line   C1         C2         P1         P2         RxDipole      PltPt     Nlevel     PltPtX     PltPtY         SP          I           Vp        Res         SD   Nstack         Mx"
        # make IP headers
        ip_columns = ""
        for ip in range(self.readings[0].Vdp[0].Vs.size):
            ip_cnt = ip + 1
            if ip_cnt < 10: 
                header_columns = header_columns + "         " + "IP0" + str(ip + 1)
            else:
                header_columns = header_columns + "         " + "IP" + str(ip + 1)

        header1 = "IP DATA FROM : " + projectName + " (TQIPdb  V2.10.432) " + date + " Mx start = " + str(start_time) + " Mx end = " + str(end_time)
        header2 = "LINE: " + line + " ARRAY:DPDP DIPOLE:" + str(a) + "  UNITS:M T=40,"
        for idx in range(len(self.window_width)):
            if idx < len(self.window_width) - 1:
                header2 = header2 + str(int(self.window_width[idx])) + ","
            else:
                header2 = header2 + str(int(self.window_width[idx]))
        # file 1
        out_file.write('%s\n' % header1)
        out_file.write('%s\n' % header2)
        out_file.write('%s\n' % header_columns)
        for rec in range(len(self.readings)):
            num_dipole = len(self.readings[rec].Vdp)
            for dp in range(num_dipole):
                SD = np.abs(self.readings[rec].Vdp[dp].Vp * self.readings[rec].Vdp[dp].Vp_err)
                vp_v = self.readings[rec].Vdp[dp].Vp / 100.
                if self.readings[rec].Vdp[dp].flagRho == "Reject":
                    out_file.write(format_dc % (self.readings[rec].Idp.Tx1East,
                                                 self.readings[rec].Idp.Tx1North,
                                                 self.readings[rec].Idp.Tx2East,
                                                 self.readings[rec].Idp.Tx2North,
                                                 self.readings[rec].Vdp[dp].Rx1East,
                                                 self.readings[rec].Vdp[dp].Rx1North,
                                                 self.readings[rec].Vdp[dp].Rx2East,
                                                 self.readings[rec].Vdp[dp].Rx2North,
                                                 self.readings[rec].Idp.Tx1x,
                                                 self.readings[rec].Idp.Tx2x,
                                                 self.readings[rec].Vdp[dp].Rx1x,
                                                 self.readings[rec].Vdp[dp].Rx2x,
                                                 self.readings[rec].Idp.Tx1y,
                                                 self.readings[rec].Idp.Tx2y,
                                                 self.readings[rec].Vdp[dp].Rx1y,
                                                 self.readings[rec].Vdp[dp].Rx2y,
                                                 self.readings[rec].Vdp[dp].getAseperationLocal(direction),
                                                 self.readings[rec].Vdp[dp].getDepthPoint2D(self.readings[rec].Idp, direction, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].getNlevel(self.readings[rec].Idp, a, direction),
                                                 self.readings[rec].Vdp[dp].getXplotpointUTM(self.readings[rec].Idp, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].getYplotpointUTM(self.readings[rec].Idp, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].Sp,
                                                 self.readings[rec].Vdp[dp].In * 1e-3,
                                                 self.readings[rec].Vdp[dp].Vp,
                                                 '*',
                                                 '*',
                                                 float(self.readings[rec].Vdp[dp].Stack),
                                                 '*'))
                    for win in range(self.readings[rec].Vdp[dp].Vs.size):
                        if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                            out_file.write('%13s' % ('*'))
                        else:
                            out_file.write('%13s\n' % ('*'))
                elif self.readings[rec].Vdp[dp].flagMx == "Reject":
                    out_file.write(format_flag % (self.readings[rec].Idp.Tx1East,
                                                 self.readings[rec].Idp.Tx1North,
                                                 self.readings[rec].Idp.Tx2East,
                                                 self.readings[rec].Idp.Tx2North,
                                                 self.readings[rec].Vdp[dp].Rx1East,
                                                 self.readings[rec].Vdp[dp].Rx1North,
                                                 self.readings[rec].Vdp[dp].Rx2East,
                                                 self.readings[rec].Vdp[dp].Rx2North,
                                                 self.readings[rec].Idp.Tx1x,
                                                 self.readings[rec].Idp.Tx2x,
                                                 self.readings[rec].Vdp[dp].Rx1x,
                                                 self.readings[rec].Vdp[dp].Rx2x,
                                                 self.readings[rec].Idp.Tx1y,
                                                 self.readings[rec].Idp.Tx2y,
                                                 self.readings[rec].Vdp[dp].Rx1y,
                                                 self.readings[rec].Vdp[dp].Rx2y,
                                                 self.readings[rec].Vdp[dp].getAseperationLocal(direction),
                                                 self.readings[rec].Vdp[dp].getDepthPoint2D(self.readings[rec].Idp, direction, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].getNlevel(self.readings[rec].Idp, a, direction),
                                                 self.readings[rec].Vdp[dp].getXplotpointUTM(self.readings[rec].Idp, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].getYplotpointUTM(self.readings[rec].Idp, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].Sp,
                                                 self.readings[rec].Vdp[dp].In * 1e-3,
                                                 self.readings[rec].Vdp[dp].Vp,
                                                 self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                                 SD,
                                                 float(self.readings[rec].Vdp[dp].Stack),
                                                 '*'))
                    for win in range(self.readings[rec].Vdp[dp].Vs.size):
                        if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                            out_file.write('%13s' % ('*'))
                        else:
                            out_file.write('%13s\n' % ('*'))
                else:
                    out_file.write(format_jup % (self.readings[rec].Idp.Tx1East,
                                                 self.readings[rec].Idp.Tx1North,
                                                 self.readings[rec].Idp.Tx2East,
                                                 self.readings[rec].Idp.Tx2North,
                                                 self.readings[rec].Vdp[dp].Rx1East,
                                                 self.readings[rec].Vdp[dp].Rx1North,
                                                 self.readings[rec].Vdp[dp].Rx2East,
                                                 self.readings[rec].Vdp[dp].Rx2North,
                                                 self.readings[rec].Idp.Tx1x,
                                                 self.readings[rec].Idp.Tx2x,
                                                 self.readings[rec].Vdp[dp].Rx1x,
                                                 self.readings[rec].Vdp[dp].Rx2x,
                                                 self.readings[rec].Idp.Tx1y,
                                                 self.readings[rec].Idp.Tx2y,
                                                 self.readings[rec].Vdp[dp].Rx1y,
                                                 self.readings[rec].Vdp[dp].Rx2y,
                                                 self.readings[rec].Vdp[dp].getAseperationLocal(direction),
                                                 self.readings[rec].Vdp[dp].getDepthPoint2D(self.readings[rec].Idp, direction, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].getNlevel(self.readings[rec].Idp, a, direction),
                                                 self.readings[rec].Vdp[dp].getXplotpointUTM(self.readings[rec].Idp, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].getYplotpointUTM(self.readings[rec].Idp, dipole_dipole=True),
                                                 self.readings[rec].Vdp[dp].Sp,
                                                 self.readings[rec].Vdp[dp].In * 1e-3,
                                                 self.readings[rec].Vdp[dp].Vp,
                                                 self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                                 SD,
                                                 float(self.readings[rec].Vdp[dp].Stack),
                                                 self.readings[rec].Vdp[dp].Mx))
                    for win in range(self.readings[rec].Vdp[dp].Vs.size):
                        if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                            out_file.write('%13.4f' % (self.readings[rec].Vdp[dp].Vs[win] / vp_v))
                        else:
                            out_file.write('%13.4f\n' % (self.readings[rec].Vdp[dp].Vs[win] / vp_v))                        

    def writeColeColeSEDat(self, outname, start_time, end_time, residuals=False):
        """
        Writes Cole-Cole data from stretched exponential to file for DiasQC
        """
        out_file = open(outname, "w+")
        if residuals is False:
            format_jup = '%6s %8s %8s %8i %13s %12.3f %12.3f %12.3f %12.0f %12.0f %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14.1f %10.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12s %14.3f %14.3f %11.0f %13s %8.3f %8.3f %8s %10.3e %10.3e %10.2f '
        else:
            format_jup = '%6s %8s %8s %8i %13s %12.3f %12.3f %12.3f %12.0f %12.0f %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14.1f %10.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12s %14.3f %14.3f %11.0f %13s %8.3f %8.3f %8s %10.3e %10.3e %10.2f %10.3e %10.3e '

        # write the headers
        for i in range(len(self.headers)):
            if i == 2:
                broken = self.headers[2].split()
                self.headers[i] = broken[0] + " " + broken[1] + " - " + str(int(start_time)) + ":" + str(int(end_time)) + "(ms)\n"
                out_file.write('%s' % self.headers[i])
            elif i < (len(self.headers) - 1):
                out_file.write('%s' % self.headers[i])
            else:
                if 'Xcorr' not in self.headers[i]:
                    broken2 = self.headers[i].split("Contact")
                    self.headers[i] = broken2[0] + "        Contact      Xcorr" + broken2[1]
                # fix headers
                hdr_temp = self.headers[i].split("Vs01")
                if residuals is False:
                    if 'Tau' in hdr_temp[0]:
                        self.headers[i] = hdr_temp[0] + "         Vs01" + hdr_temp[1]
                    else:
                        self.headers[i] = hdr_temp[0] + "   C          Tau          M         Vs01" + hdr_temp[1]
                else:
                    if 'rdc' in hdr_temp[0]:
                        self.headers[i] = hdr_temp[0] + "         Vs01" + hdr_temp[1]
                    elif 'Tau' in hdr_temp[0]:
                        self.headers[i] = hdr_temp[0] + "         rdc         rip         Vs01" + hdr_temp[1]
                    else:
                        self.headers[i] = hdr_temp[0] + "   C          Tau          M         rdc         rip         Vs01" + hdr_temp[1]
                out_file.write('%s' % self.headers[i])

        for rec in range(len(self.readings)):
            num_dipole = len(self.readings[rec].Vdp)
            for dp in range(num_dipole):
                if self.readings[rec].Vdp[dp].c is None:
                    cole_con = -99.9
                    cole_tau = -99.9
                    cole_eta = -99.9
                else:
                    cole_con = self.readings[rec].Vdp[dp].c
                    cole_tau = self.readings[rec].Vdp[dp].tau
                    cole_eta = self.readings[rec].Vdp[dp].eta
                if np.abs(self.readings[rec].Vdp[dp].Mx_err) > 999.0:
                    self.readings[rec].Vdp[dp].Mx_err = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                if np.abs(self.readings[rec].Vdp[dp].Mx) > 999.0:
                    self.readings[rec].Vdp[dp].Mx = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                if residuals is False:
                    out_file.write(format_jup % (str(self.readings[rec].Vdp[dp].reading),
                                                 str(self.readings[rec].Vdp[dp].dipole),
                                                 str(self.readings[rec].Vdp[dp].status),
                                                 int(self.readings[rec].Vdp[dp].getAseperation()),
                                                 self.readings[rec].Idp.Tx1File,
                                                 self.readings[rec].Idp.Tx1East,
                                                 self.readings[rec].Idp.Tx1North,
                                                 self.readings[rec].Idp.Tx1Elev,
                                                 self.readings[rec].Idp.Tx1x,
                                                 self.readings[rec].Idp.Tx1y,
                                                 self.readings[rec].Idp.Tx2East,
                                                 self.readings[rec].Idp.Tx2North,
                                                 self.readings[rec].Idp.Tx2Elev,
                                                 self.readings[rec].Idp.Tx2x,
                                                 self.readings[rec].Idp.Tx2y,
                                                 self.readings[rec].Vdp[dp].Rx1File,
                                                 self.readings[rec].Vdp[dp].Rx1East,
                                                 self.readings[rec].Vdp[dp].Rx1North,
                                                 self.readings[rec].Vdp[dp].Rx1Elev,
                                                 self.readings[rec].Vdp[dp].Rx1x,
                                                 self.readings[rec].Vdp[dp].Rx1y,
                                                 self.readings[rec].Vdp[dp].Rx2File,
                                                 self.readings[rec].Vdp[dp].Rx2East,
                                                 self.readings[rec].Vdp[dp].Rx2North,
                                                 self.readings[rec].Vdp[dp].Rx2Elev,
                                                 self.readings[rec].Vdp[dp].Rx2x,
                                                 self.readings[rec].Vdp[dp].Rx2y,
                                                 0.0,
                                                 self.readings[rec].Vdp[dp].xcorr,
                                                 self.readings[rec].Vdp[dp].Sp,
                                                 self.readings[rec].Vdp[dp].Vp,
                                                 self.readings[rec].Vdp[dp].Vp_err,
                                                 self.readings[rec].Vdp[dp].In,
                                                 self.readings[rec].Vdp[dp].In_err,
                                                 self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                                 self.readings[rec].Vdp[dp].flagRho,
                                                 self.readings[rec].Vdp[dp].calcGeoFactor(self.readings[rec].Idp),
                                                 self.readings[rec].Vdp[dp].coupling,
                                                 self.readings[rec].Vdp[dp].Stack,
                                                 self.readings[rec].Vdp[dp].TimeBase,
                                                 self.readings[rec].Vdp[dp].Mx,
                                                 self.readings[rec].Vdp[dp].Mx_err,
                                                 self.readings[rec].Vdp[dp].flagMx,
                                                 cole_con,
                                                 cole_tau,
                                                 cole_eta))
                else:
                    out_file.write(format_jup % (str(self.readings[rec].Vdp[dp].reading),
                                                 str(self.readings[rec].Vdp[dp].dipole),
                                                 str(self.readings[rec].Vdp[dp].status),
                                                 int(self.readings[rec].Vdp[dp].getAseperation()),
                                                 self.readings[rec].Idp.Tx1File,
                                                 self.readings[rec].Idp.Tx1East,
                                                 self.readings[rec].Idp.Tx1North,
                                                 self.readings[rec].Idp.Tx1Elev,
                                                 self.readings[rec].Idp.Tx1x,
                                                 self.readings[rec].Idp.Tx1y,
                                                 self.readings[rec].Idp.Tx2East,
                                                 self.readings[rec].Idp.Tx2North,
                                                 self.readings[rec].Idp.Tx2Elev,
                                                 self.readings[rec].Idp.Tx2x,
                                                 self.readings[rec].Idp.Tx2y,
                                                 self.readings[rec].Vdp[dp].Rx1File,
                                                 self.readings[rec].Vdp[dp].Rx1East,
                                                 self.readings[rec].Vdp[dp].Rx1North,
                                                 self.readings[rec].Vdp[dp].Rx1Elev,
                                                 self.readings[rec].Vdp[dp].Rx1x,
                                                 self.readings[rec].Vdp[dp].Rx1y,
                                                 self.readings[rec].Vdp[dp].Rx2File,
                                                 self.readings[rec].Vdp[dp].Rx2East,
                                                 self.readings[rec].Vdp[dp].Rx2North,
                                                 self.readings[rec].Vdp[dp].Rx2Elev,
                                                 self.readings[rec].Vdp[dp].Rx2x,
                                                 self.readings[rec].Vdp[dp].Rx2y,
                                                 0.0,
                                                 self.readings[rec].Vdp[dp].xcorr,
                                                 self.readings[rec].Vdp[dp].Sp,
                                                 self.readings[rec].Vdp[dp].Vp,
                                                 self.readings[rec].Vdp[dp].Vp_err,
                                                 self.readings[rec].Vdp[dp].In,
                                                 self.readings[rec].Vdp[dp].In_err,
                                                 self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                                 self.readings[rec].Vdp[dp].flagRho,
                                                 self.readings[rec].Vdp[dp].calcGeoFactor(self.readings[rec].Idp),
                                                 self.readings[rec].Vdp[dp].coupling,
                                                 self.readings[rec].Vdp[dp].Stack,
                                                 self.readings[rec].Vdp[dp].TimeBase,
                                                 self.readings[rec].Vdp[dp].Mx,
                                                 self.readings[rec].Vdp[dp].Mx_err,
                                                 self.readings[rec].Vdp[dp].flagMx,
                                                 cole_con,
                                                 cole_tau,
                                                 cole_eta,
                                                 self.readings[rec].Vdp[dp].residual_dc,
                                                 self.readings[rec].Vdp[dp].residual_ip))
                for win in range(self.readings[rec].Vdp[dp].Vs.size):
                    if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                        out_file.write('%12.4f' % (self.readings[rec].Vdp[dp].Vs[win]))
                    else:
                        out_file.write('%12.4f\n' % (self.readings[rec].Vdp[dp].Vs[win]))

    def writeColeColeDat(self, eta, tau, c,
                         error, indicies, outname, start_time, end_time):
        """
        Writes Cole-Cole data from stretched exponential to file for DiasQC
        """
        for idx in range(indicies.shape[0]):
            rdg = indicies[idx, 0]
            dp = indicies[idx, 1]
            self.readings[rdg].Vdp[dp].c = c[idx]
            self.readings[rdg].Vdp[dp].eta = eta[idx]
            self.readings[rdg].Vdp[dp].tau = tau[idx]
            if error[idx] == -99.9:
                self.readings[rdg].Vdp[dp].Mx_err = error[idx]
                self.readings[rdg].Vdp[dp].flagMx = "Reject"
            else:
                self.readings[rdg].Vdp[dp].Mx_err = error[idx]
                self.readings[rdg].Vdp[dp].flagMx = "Accept"

        out_file = open(outname, "w+")
        format_jup = '%6s %8s %8s %8i %13s %12.3f %12.3f %12.3f %12.0f %12.0f %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14s %12.3f %12.3f %12.3f %12.0f %12.0f %14.1f %10.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12s %14.3f %14.3f %11.0f %13s %8.3f %8.3f %8s %12.3e %12.3e %10.2f '
        # write the headers
        for i in range(len(self.headers)):
            if i == 2:
                broken = self.headers[2].split()
                self.headers[i] = broken[0] + " " + broken[1] + " - " + str(int(start_time)) + ":" + str(int(end_time)) + "(ms)\n"
                out_file.write('%s' % self.headers[i])
            elif i < (len(self.headers) - 1):
                out_file.write('%s' % self.headers[i])
            else:
                # fix headers
                hdr_temp = self.headers[i].split("Vs01")
                self.headers[i] = hdr_temp[0] + "   C          Tau          M         Vs01" + hdr_temp[1]
                out_file.write('%s' % self.headers[i])

        for rec in range(len(self.readings)):
            num_dipole = len(self.readings[rec].Vdp)
            for dp in range(num_dipole):
                if self.readings[rec].Vdp[dp].c is None:
                    cole_con = -99.9
                    cole_tau = -99.9
                    cole_eta = -99.9
                else:
                    cole_con = self.readings[rec].Vdp[dp].c
                    cole_tau = self.readings[rec].Vdp[dp].tau
                    cole_eta = self.readings[rec].Vdp[dp].eta
                if np.abs(self.readings[rec].Vdp[dp].Mx_err) > 999.0:
                    self.readings[rec].Vdp[dp].Mx_err = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                if np.abs(self.readings[rec].Vdp[dp].Mx) > 999.0:
                    self.readings[rec].Vdp[dp].Mx = -99.9
                    self.readings[rec].Vdp[dp].flagMx = "Reject"
                out_file.write(format_jup % (str(self.readings[rec].Vdp[dp].reading),
                                             str(self.readings[rec].Vdp[dp].dipole),
                                             str(self.readings[rec].Vdp[dp].status),
                                             int(self.readings[rec].Vdp[dp].getAseperation()),
                                             self.readings[rec].Idp.Tx1File,
                                             self.readings[rec].Idp.Tx1East,
                                             self.readings[rec].Idp.Tx1North,
                                             self.readings[rec].Idp.Tx1Elev,
                                             self.readings[rec].Idp.Tx1x,
                                             self.readings[rec].Idp.Tx1y,
                                             self.readings[rec].Idp.Tx2East,
                                             self.readings[rec].Idp.Tx2North,
                                             self.readings[rec].Idp.Tx2Elev,
                                             self.readings[rec].Idp.Tx2x,
                                             self.readings[rec].Idp.Tx2y,
                                             self.readings[rec].Vdp[dp].Rx1File,
                                             self.readings[rec].Vdp[dp].Rx1East,
                                             self.readings[rec].Vdp[dp].Rx1North,
                                             self.readings[rec].Vdp[dp].Rx1Elev,
                                             self.readings[rec].Vdp[dp].Rx1x,
                                             self.readings[rec].Vdp[dp].Rx1y,
                                             self.readings[rec].Vdp[dp].Rx2File,
                                             self.readings[rec].Vdp[dp].Rx2East,
                                             self.readings[rec].Vdp[dp].Rx2North,
                                             self.readings[rec].Vdp[dp].Rx2Elev,
                                             self.readings[rec].Vdp[dp].Rx2x,
                                             self.readings[rec].Vdp[dp].Rx2y,
                                             0.0,
                                             self.readings[rec].Vdp[dp].xcorr,
                                             self.readings[rec].Vdp[dp].Sp,
                                             self.readings[rec].Vdp[dp].Vp,
                                             self.readings[rec].Vdp[dp].Vp_err,
                                             self.readings[rec].Vdp[dp].In,
                                             self.readings[rec].Vdp[dp].In_err,
                                             self.readings[rec].Vdp[dp].calcRho(self.readings[rec].Idp),
                                             self.readings[rec].Vdp[dp].flagRho,
                                             self.readings[rec].Vdp[dp].calcGeoFactor(self.readings[rec].Idp),
                                             self.readings[rec].Vdp[dp].coupling,
                                             self.readings[rec].Vdp[dp].Stack,
                                             self.readings[rec].Vdp[dp].TimeBase,
                                             self.readings[rec].Vdp[dp].Mx,
                                             self.readings[rec].Vdp[dp].Mx_err,
                                             self.readings[rec].Vdp[dp].flagMx,
                                             cole_con,
                                             cole_tau,
                                             cole_eta))
                for win in range(self.readings[rec].Vdp[dp].Vs.size):
                    if win < (self.readings[rec].Vdp[dp].Vs.size - 1):
                        out_file.write('%12.4f' % (self.readings[rec].Vdp[dp].Vs[win]))
                    else:
                        out_file.write('%12.4f\n' % (self.readings[rec].Vdp[dp].Vs[win]))

    def getApparentResistivity(self, reject="app_rho", calc=True):
        """
        Exports all the apparent resistivity data

        Output:
        numpy array [value]

        """
        resistivity_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        if calc:
                            rho_a = self.readings[k].Vdp[j].calcRho(self.readings[k].Idp)
                            resistivity_list.append(rho_a)
                        else:
                            rho_a = self.readings[k].Vdp[j].Rho
                            resistivity_list.append(rho_a)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        if calc:
                            rho_a = self.readings[k].Vdp[j].calcRho(self.readings[k].Idp)
                            resistivity_list.append(rho_a)
                        else:
                            rho_a = self.readings[k].Vdp[j].Rho
                            resistivity_list.append(rho_a)

        return np.vstack(resistivity_list)

    def getResidualsXyz(self, data_type='dc', file_output=None):
        """
        Exports all the frequency component data

        Output:
        numpy array [value]

        """
        delta_list = []
        x_list = []
        y_list = []
        z_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if data_type == 'dc':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        delta_list.append(self.readings[k].Vdp[j].residual_dc)
                        if file_output is not None:
                            x_list.append(self.readings[k].Vdp[j].getXplotpointUTM(self.readings[k].Idp))
                            y_list.append(self.readings[k].Vdp[j].getYplotpointUTM(self.readings[k].Idp))
                            z_list.append(self.readings[k].Vdp[j].getZplotpoint(self.readings[k].Idp))
                elif data_type == 'ip':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        delta_list.append(self.readings[k].Vdp[j].residual_ip)
                        if file_output is not None:
                            x_list.append(self.readings[k].Vdp[j].getXplotpointUTM(self.readings[k].Idp))
                            y_list.append(self.readings[k].Vdp[j].getYplotpointUTM(self.readings[k].Idp))
                            z_list.append(self.readings[k].Vdp[j].getZplotpoint(self.readings[k].Idp))
                else:
                    print('[WARNING] data type requested is invalid')
        if file_output is not None:
            out_file = open(file_output, "w+")
            for idx in range(len(delta_list)):
                if x_list[idx] < 0 or y_list[idx] < 0:
                    pass
                else:
                    out_file.write('%0.2f,%0.2f,%0.2f,%0.4e\n' % (x_list[idx], y_list[idx], z_list[idx], delta_list[idx]))
            out_file.close()
        return np.vstack(delta_list)

    def getFrequencyComponent(self, reject="app_rho"):
        """
        Exports all the frequency component data

        Output:
        numpy array [value]

        """
        c_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        c_list.append(self.readings[k].Vdp[j].c)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        c_list.append(self.readings[k].Vdp[j].c)

        return np.vstack(c_list)

    def getTimeConstants(self, reject="app_rho"):
        """
        Exports all the time constant data

        Output:
        numpy array [value]

        """
        tau_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        tau_list.append(self.readings[k].Vdp[j].tau)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        tau_list.append(self.readings[k].Vdp[j].tau)

        return np.vstack(tau_list)

    def getColeColeMx(self, reject="Rho"):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        eta_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        eta_list.append(self.readings[k].Vdp[j].eta)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        eta_list.append(self.readings[k].Vdp[j].eta)

        return np.vstack(eta_list)

    def getVsErrors(self, reject=None):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        error_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        error_list.append(self.readings[k].Vdp[j].Mx_err)
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        error_list.append(self.readings[k].Vdp[j].Mx_err)
                else:
                    error_list.append(self.readings[k].Vdp[j].Mx_err)

        return np.vstack(error_list)

    def getDCvoltages(self, reject=None):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        vp_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        vp_list.append(self.readings[k].Vdp[j].Vp / 1e3)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                       vp_list.append(self.readings[k].Vdp[j].Vp / 1e3)
                else:
                    vp_list.append(self.readings[k].Vdp[j].Vp / 1e3)

        return np.vstack(vp_list)

    def getVpDivInErrors(self, reject=None):
        """
        Exports all the eta data

        Output:
        numpy array [value]

        """
        error_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        err = (self.readings[k].Vdp[j].Vp_err + self.readings[k].Vdp[j].In_err) / 100.
                        error_list.append(err)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        err = (self.readings[k].Vdp[j].Vp_err + self.readings[k].Vdp[j].In_err) / 100.
                        error_list.append(err)
                else:
                    err = (self.readings[k].Vdp[j].Vp_err + self.readings[k].Vdp[j].In_err) / 100.
                    error_list.append(err)

        return np.vstack(error_list)

    def getActiveIndicies(self, reject=None):
        """
        exports an array containing index of each available dipole

        Output:
        list []
        """
        index_list = []
        num_rdg = len(self.readings)
        if reject is not None:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject is 'app_rho':
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            index = np.array([k, j])
                            index_list.append(index)
                    elif reject is 'app_mx':
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            index = np.array([k, j])
                            index_list.append(index)
        else:
            print("Please choose which indicies you are looking for")

        return np.vstack(index_list)

    def getPercentageNegativeAppRho(self, reject=None):
        """
        calculates % of -ve app rho values

        Output:
        numpy array [value]

        """
        cnt_neg_app_rho = 0
        total_dipoles = 0
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        total_dipoles += 1
                        if self.readings[k].Vdp[j].Rho < 0:
                            cnt_neg_app_rho += 1
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        total_dipoles += 1
                        if self.readings[k].Vdp[j].Rho < 0:
                            cnt_neg_app_rho += 1
                else:
                    total_dipoles += 1
                    if self.readings[k].Vdp[j].Rho < 0:
                            cnt_neg_app_rho += 1
        percent_negative = (cnt_neg_app_rho / total_dipoles) * 100
        print("{1} - Percentage of -ve apparent Resistivities: {0} %".format(percent_negative, reject))

    def getPercentageNegativeAppMx(self, reject=None):
        """
        calculates % of -ve app rho values

        Output:
        numpy array [value]

        """
        cnt_neg_app_mx = 0
        total_dipoles = 0
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        total_dipoles += 1
                        if self.readings[k].Vdp[j].Mx < 0:
                            cnt_neg_app_mx += 1
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        total_dipoles += 1
                        if self.readings[k].Vdp[j].Mx < 0:
                            cnt_neg_app_mx += 1
                else:
                    total_dipoles += 1
                    if self.readings[k].Vdp[j].Mx < 0:
                            cnt_neg_app_mx += 1
        percent_negative = (cnt_neg_app_mx / total_dipoles) * 100
        print("{1} - Percentage of -ve apparent Mx: {0} %".format(percent_negative, reject))

    def getPercentageLowVoltages(self, threshold=1.0, reject=None):
        """
        calculates % of -ve app rho values

        Output:
        numpy array [value]

        """
        cnt_low_vp = 0
        total_dipoles = 0
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject is "app_rho":
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        total_dipoles += 1
                        if np.abs(self.readings[k].Vdp[j].Vp) < threshold:
                            cnt_low_vp += 1
                elif reject is "app_mx":
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        total_dipoles += 1
                        if np.abs(self.readings[k].Vdp[j].Vp) < threshold:
                            cnt_low_vp += 1
                else:
                    total_dipoles += 1
                    if np.abs(self.readings[k].Vdp[j].Vp) < threshold:
                            cnt_low_vp += 1
        percent_low = (cnt_low_vp / total_dipoles) * 100
        print("{2} - Percentage of Voltages under {0} mV are: {1} %".format(threshold, percent_low, reject))

    def getGeometricFactor(self, reject=None):
        """
        Exports all the geometry factor data

        Output:
        numpy array [value]

        """
        k_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_rho':
                    if (self.readings[k].Vdp[j].flagRho == "Accept"):
                        k_a = self.readings[k].Vdp[j].calcGeoFactor(self.readings[k].Idp)
                        k_list.append(k_a)
                elif reject == 'app_mx':
                    if (self.readings[k].Vdp[j].flagMx == "Accept"):
                        k_a = self.readings[k].Vdp[j].calcGeoFactor(self.readings[k].Idp)
                        k_list.append(k_a)
                else:
                    k_a = self.readings[k].Vdp[j].calcGeoFactor(self.readings[k].Idp)
                    k_list.append(k_a)
        return np.asarray(k_list)

    def getVoltages(self, reject=None):
        """
        Exports all the geometry factor data

        Output:
        numpy array [value]

        """
        vp_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        Vp = self.readings[k].Vdp[j].Vp
                        vp_list.append(Vp)
                        if Vp == 0:
                            print("[WARNING] Vp of ZERO at rdg: {0} & dp: {1}".format(self.readings[k].Vdp[j].reading, self.readings[k].Vdp[j].dipole))
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        Vp = self.readings[k].Vdp[j].Vp
                        vp_list.append(Vp)
                        if Vp == 0:
                            print("[WARNING] Vp of ZERO at rdg: {0} & dp: {1}".format(self.readings[k].Vdp[j].reading, self.readings[k].Vdp[j].dipole)) 
                else:
                    Vp = self.readings[k].Vdp[j].Vp
                    vp_list.append(Vp)
                    if Vp == 0:
                        print("[WARNING] Vp of ZERO at rdg: {0} & dp: {1}".format(self.readings[k].Vdp[j].reading, self.readings[k].Vdp[j].dipole)) 

        return np.asarray(vp_list)

    def reCalcApparentChageability(self, start, stop, reject=None):
        """
        Exports all the apparent chargeability data

        Output:
        numpy array [value]

        """
        chargeability_list = []
        num_rdg = len(self.readings)
        window_widths = self.window_width
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        mx_a = self.readings[k].Vdp[j].Mx
                        if not np.isinf(mx_a):
                            mx_a = self.readings[k].Vdp[j].calcMx(window_widths, start, stop)
                            chargeability_list.append(mx_a)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        mx_a = self.readings[k].Vdp[j].Mx
                        if not np.isinf(mx_a):
                            mx_a = self.readings[k].Vdp[j].calcMx(window_widths, start, stop)
                            chargeability_list.append(mx_a)
                else:
                    mx_a = self.readings[k].Vdp[j].Mx
                    if not np.isinf(mx_a):
                        mx_a = self.readings[k].Vdp[j].calcMx(window_widths, start, stop)
                        chargeability_list.append(mx_a)

        return np.asarray(chargeability_list)

    def getApparentChageability(self, reject=None):
        """
        Exports all the apparent chargeability data

        Output:
        numpy array [value]

        """
        chargeability_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        mx_a = self.readings[k].Vdp[j].Mx
                        if np.isinf(mx_a):
                            print("[WARNING] found inf Mx at rdg: {0} & dp: {1}".format(self.readings[k].Vdp[j].reading, self.readings[k].Vdp[j].dipole))
                        else:
                            chargeability_list.append(mx_a)
                elif reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        mx_a = self.readings[k].Vdp[j].Mx
                        if np.isinf(mx_a):
                            print("[WARNING] found inf Mx at rdg: {0} & dp: {1}".format(self.readings[k].Vdp[j].reading, self.readings[k].Vdp[j].dipole))
                        else:
                            chargeability_list.append(mx_a)
                else:
                    mx_a = self.readings[k].Vdp[j].Mx
                    if np.isinf(mx_a):
                        print("[WARNING] found inf Mx at rdg: {0} & dp{1}: ".format(self.readings[k].Vdp[j].reading, self.readings[k].Vdp[j].dipole))
                    else:
                        chargeability_list.append(mx_a)

        return np.asarray(chargeability_list)

    def getNullFactors(self, reject=None):
        """
        Exports all the coupling factor data

        Output:
        numpy array [value]

        """
        null_list = []
        num_rdg = len(self.readings)
        for rdg in range(num_rdg):
            num_dipole = len(self.readings[rdg].Vdp)
            for dp in range(num_dipole):
                if reject == 'app_rho':
                    if self.readings[rdg].Vdp[dp].flagRho == "Accept":
                        null_list.append(self.readings[rdg].Vdp[dp].coupling)
                elif reject == 'app_mx':
                    if self.readings[rdg].Vdp[dp].flagMx == "Accept":
                        null_list.append(self.readings[rdg].Vdp[dp].coupling)
                else:
                        null_list.append(self.readings[rdg].Vdp[dp].coupling)

        return np.asarray(null_list)

    def getCurrents(self, reject=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        in_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if reject == 'app_rho':
                    if self.readings[k].Vdp[j].flagRho == "Accept":
                        In = self.readings[k].Vdp[j].In
                        in_list.append(In)
                elif reject == 'app_mx':
                    if self.readings[k].Vdp[j].flagMx == "Accept":
                        In = self.readings[k].Vdp[j].In
                        in_list.append(In)
                else:
                    In = self.readings[k].Vdp[j].In
                    in_list.append(In)

        return np.asarray(in_list)

    def calculateDpredMinusData(self, dpred=None, data_type=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        inds = self.getActiveIndicies(reject='app_rho')
        if data_type.lower() == 'dc':
            inds = self.getActiveIndicies(reject='app_rho')
        else:
            inds = self.getActiveIndicies(reject='app_mx')
        print('size of data: {0}, inds: {1}'.format(dpred.size, inds.shape))
        for idx in range(dpred.size):
            if data_type.lower() == 'ip':
                Vp = self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].Vp * 1e-3
                In = self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].In * 1e-3
                voltage = (self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].Mx * 1e-3) * (Vp / In)
                residual = dpred[idx] - voltage
                self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].residual_ip = residual
            if data_type.lower() == 'dc':
                resistance = self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].Vp / self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].In
                residual = dpred[idx] - resistance
                self.readings[inds[idx, 0]].Vdp[inds[idx, 1]].residual_dc = residual

    def getAspacings(self, local=False, direction=None, reject=None):
        a_list = []
        num_rdg = len(self.readings)
        if local:
            if direction is not None:
                for k in range(num_rdg):
                    num_dipole = len(self.readings[k].Vdp)
                    for j in range(num_dipole):
                        if reject == 'app_rho':
                            if self.readings[k].Vdp[j].flagRho == "Accept":
                                sep = np.abs(self.readings[k].Vdp[j].getAseperationLocal(direction=direction))
                                a_list.append(sep)
                        if reject == 'app_mx':
                            if self.readings[k].Vdp[j].flagMx == "Accept":
                                sep = np.abs(self.readings[k].Vdp[j].getAseperationLocal(direction=direction))
                                a_list.append(sep)
                        else:
                            sep = np.abs(self.readings[k].Vdp[j].getAseperationLocal(direction=direction))
                            a_list.append(sep)
            else:
                print("[OUTPUT]  Please define line direction for local option!!!!!!!!!")
        else:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject == 'app_rho':
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                            a_list.append(sep)
                    elif reject == 'app_mx':
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                            a_list.append(sep)
                    else:
                        sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                        a_list.append(sep)
        return a_list

    def getIndicesofAspacing(self, spacing_min=None, spacing_max=None, reject=None):
        a_list = []
        num_rdg = len(self.readings)
        if spacing_min is not None:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject == 'app_rho':
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                            if spacing_min <= sep <= spacing_max:
                                a_list.append(True)
                            else:
                                a_list.append(False)
                    elif reject == 'app_mx':
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                            if spacing_min <= sep <= spacing_max:
                                a_list.append(True)
                            else:
                                a_list.append(False)
                    else:
                        sep = np.abs(self.readings[k].Vdp[j].getAseperation())
                        a_list.append(sep)
        else:
            print('[WARNING] Pick an a spacing! this is what the function is for!!!')
        return a_list

    def getIndicesofCurrentRange(self, cutoff=None, reject=None):
        a_list = []
        num_rdg = len(self.readings)
        if cutoff is not None:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject == 'app_rho':
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            if self.readings[k].Vdp[j].In <= cutoff:
                                a_list.append(True)
                            else:
                                a_list.append(False)
                    elif reject == 'app_mx':
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            if self.readings[k].Vdp[j].In <= cutoff:
                                a_list.append(True)
                            else:
                                a_list.append(False)
                    else:
                        if self.readings[k].Vdp[j].In <= cutoff:
                            a_list.append(True)
                        else:
                            a_list.append(False)
        else:
            print('[WARNING] Pick an a spacing! this is what the function is for!!!')
        return a_list

    def getIndicesofTxRxspacing(self, spacing_min=None, spacing_max=None, reject=None):
        a_list = []
        num_rdg = len(self.readings)
        if spacing_min is not None:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    if reject == 'app_rho':
                        if self.readings[k].Vdp[j].flagRho == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getTxRxSeperation(self.readings[k].Idp))
                            if spacing_min <= sep <= spacing_max:
                                a_list.append(True)
                            else:
                                a_list.append(False)
                    elif reject == 'app_mx':
                        if self.readings[k].Vdp[j].flagMx == "Accept":
                            sep = np.abs(self.readings[k].Vdp[j].getTxRxSeperation(self.readings[k].Idp))
                            if spacing_min <= sep <= spacing_max:
                                a_list.append(True)
                            else:
                                a_list.append(False)
                    else:
                        sep = np.abs(self.readings[k].Vdp[j].getTxRxSeperation(self.readings[k].Idp))
                        a_list.append(sep)
        else:
            print('[WARNING] Pick an a spacing! this is what the function is for!!!')
        return a_list

    def getSources(self, dipole=False, local=False, line_dir=None):
        """
        Exports all the tx locations to a numpy array

        Output:
        numpy array [x0,y0,z0,x1,y1,z1] (c1,c2)

        """
        src_list = []
        num_rdg = len(self.readings)
        if dipole == False:
            for k in range(num_rdg):
                if local == False:
                    tx = np.array([self.readings[k].Idp.Tx1East,
                                   self.readings[k].Idp.Tx1North,
                                   self.readings[k].Idp.Tx1Elev,
                                   self.readings[k].Idp.Tx2East,
                                   self.readings[k].Idp.Tx2North,
                                   self.readings[k].Idp.Tx2Elev])
                    src_list.append(tx)
                else:
                    if line_dir == 'ns':
                        tx = np.array([self.readings[k].Idp.Tx1y,
                               self.readings[k].Idp.Tx1Elev,
                               self.readings[k].Idp.Tx2y,
                               self.readings[k].Idp.Tx2Elev])
                    elif line_dir == 'ew':
                        tx = np.array([self.readings[k].Idp.Tx1x,
                               self.readings[k].Idp.Tx1Elev,
                               self.readings[k].Idp.Tx2x,
                               self.readings[k].Idp.Tx2Elev])
                    else:
                        print('[ERROR] Line dir must be ns - ew!!')
                        tx = np.array([self.readings[k].Idp.Tx1x,
                               self.readings[k].Idp.Tx1y,
                               self.readings[k].Idp.Tx1Elev,
                               self.readings[k].Idp.Tx2x,
                               self.readings[k].Idp.Tx2y,
                               self.readings[k].Idp.Tx2Elev])
                    src_list.append(tx)
        else:
            for k in range(num_rdg):
                num_dipole = len(self.readings[k].Vdp)
                for j in range(num_dipole):
                    print(line_dir, local)
                    if local == False:
                        tx = np.array([self.readings[k].Idp.Tx1East,
                                       self.readings[k].Idp.Tx1North,
                                       self.readings[k].Idp.Tx1Elev,
                                       self.readings[k].Idp.Tx2East,
                                       self.readings[k].Idp.Tx2North,
                                       self.readings[k].Idp.Tx2Elev])
                        src_list.append(tx)
                    else:
                        if line_dir is None:
                            tx = np.array([self.readings[k].Idp.Tx1x,
                                       self.readings[k].Idp.Tx1y,
                                       self.readings[k].Idp.Tx1Elev,
                                       self.readings[k].Idp.Tx2x,
                                       self.readings[k].Idp.Tx2y,
                                       self.readings[k].Idp.Tx2Elev])
                        else:
                            if line_dir == 'ns':
                                tx = np.array([self.readings[k].Idp.Tx1y,
                                       self.readings[k].Idp.Tx1Elev,
                                       self.readings[k].Idp.Tx2y,
                                       self.readings[k].Idp.Tx2Elev])
                            elif line_dir == 'ew':
                                tx = np.array([self.readings[k].Idp.Tx1x,
                                       self.readings[k].Idp.Tx1Elev,
                                       self.readings[k].Idp.Tx2x,
                                       self.readings[k].Idp.Tx2Elev])
                            else:
                                print('[ERROR] Line dir must be ns - ew!!')
                                tx = np.array([self.readings[k].Idp.Tx1x,
                               self.readings[k].Idp.Tx1y,
                               self.readings[k].Idp.Tx1Elev,
                               self.readings[k].Idp.Tx2x,
                               self.readings[k].Idp.Tx2y,
                               self.readings[k].Idp.Tx2Elev])
                        src_list.append(tx)

        return np.asarray(src_list)

    def getDipoles(self, local=False, line_dir=None):
        """
        Exports all the tx locations to a numpy array

        Output:
        numpy array [x0,y0,z0,x1,y1,z1] (p1,p2)

        """
        dipole_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if local == False:
                    rx = np.array([self.readings[k].Vdp[j].Rx1East,
                                   self.readings[k].Vdp[j].Rx1North,
                                   self.readings[k].Vdp[j].Rx1Elev,
                                   self.readings[k].Vdp[j].Rx2East,
                                   self.readings[k].Vdp[j].Rx2North,
                                   self.readings[k].Vdp[j].Rx2Elev])
                    dipole_list.append(rx)
                else:
                    if line_dir is None:
                        rx = np.array([self.readings[k].Vdp[j].Rx1x,
                                       self.readings[k].Vdp[j].Rx1y,
                                       self.readings[k].Vdp[j].Rx1Elev,
                                       self.readings[k].Vdp[j].Rx2x,
                                       self.readings[k].Vdp[j].Rx2y,
                                       self.readings[k].Vdp[j].Rx2Elev])
                    else:
                        if line_dir == 'ns':
                            rx = np.array([self.readings[k].Vdp[j].Rx1y,
                                   self.readings[k].Vdp[j].Rx1Elev,
                                   self.readings[k].Vdp[j].Rx2y,
                                   self.readings[k].Vdp[j].Rx2Elev])
                        elif line_dir == 'ew':
                            rx = np.array([self.readings[k].Vdp[j].Rx1x,
                                   self.readings[k].Vdp[j].Rx1Elev,
                                   self.readings[k].Vdp[j].Rx2x,
                                   self.readings[k].Vdp[j].Rx2Elev])
                        else:
                            print('[ERROR] Line dir must be ns - ew!!')
                            return
                    dipole_list.append(rx)
        number_of_dipoles = len(dipole_list)
        if line_dir is None:
            dipoles = np.zeros((number_of_dipoles, 6))
        else:
            dipoles = np.zeros((number_of_dipoles, 4))
        for rx in range(number_of_dipoles):
            dipoles[rx, :] = dipole_list[rx]
        return dipoles

    def getNumberOfDataPoints(self, active=False):
        """
        returns number of dipoles in the dataset
        """
        num_rdg = len(self.readings)
        count = 0
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                # if active:
                #     if self.readings[k].Vdp[j].flagRho
                count += 1
        return count

    def CompareDatabase2RecalcValues(self, mx_start, mx_end, reject=None):
        calc_mx = self.reCalcApparentChageability(mx_start,
                                                  mx_end, reject=reject)
        mx = self.getApparentChageability(reject=reject)
        calc_rho = self.getApparentResistivity(reject=reject)
        rho = self.getApparentResistivity(reject=reject, calc=False)
        fig = plt.figure()
        fig.suptitle('Calculated Vs. Database', fontsize=16)
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        mean_drho = str(np.mean(calc_rho - rho))
        mean_dmx = str(np.mean(calc_mx - mx))

        ax3.set_title(r'$\rho$')
        ax3.set_ylabel(r'$\Omega$' + '-m')
        ax4.set_title(r'$\Delta\rho$' + '     mean: ' + mean_drho)
        ax1.set_title(r'$\eta$')
        ax1.set_ylabel('mV/V')
        ax2.set_title(r'$\Delta\eta$' + '     mean: ' + mean_dmx)

        ax1.plot(mx, '.')
        ax1.plot(calc_mx, '.r', markersize=1)
        ax2.plot(calc_mx - mx, '.g', markersize=1)
        ax3.plot(rho, '.')
        ax3.plot(calc_rho, '.r', markersize=1)
        ax4.plot(calc_rho - rho, '.g', markersize=1)
        plt.show()

    def rejectMx(self, minv=None, maxv=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        in_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if self.readings[k].Vdp[j].flagMx == "Accept":
                    if self.readings[k].Vdp[j].Mx < minv:
                        self.readings[k].Vdp[j].flagMx = "Reject"
                    if self.readings[k].Vdp[j].Mx > maxv:
                        self.readings[k].Vdp[j].flagMx = "Reject"

        return np.asarray(in_list)

    def rejectIPresidual(self, minv=None, maxv=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if self.readings[k].Vdp[j].flagMx == "Accept":
                    if self.readings[k].Vdp[j].residual_ip < minv:
                        self.readings[k].Vdp[j].flagMx = "Reject"
                    if self.readings[k].Vdp[j].residual_ip > maxv:
                        self.readings[k].Vdp[j].flagMx = "Reject"

    def getIPresiduals(self, uncert=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        in_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if self.readings[k].Vdp[j].flagMx == "Accept":
                    in_list.append(self.readings[k].Vdp[j].residual_ip)
        residuals = np.asarray(in_list)
        # normalise the residuals if data is supplied
        if uncert is not None:
            residuals = residuals / uncert
        return residuals

    def getDCresiduals(self, uncert=None):
        """
        Exports all the In data

        Output:
        numpy array [value]

        """
        in_list = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            for j in range(num_dipole):
                if self.readings[k].Vdp[j].flagRho == "Accept":
                    in_list.append(self.readings[k].Vdp[j].residual_dc)
        residuals = np.asarray(in_list)
        # normalise the residuals if data is supplied
        if uncert is not None:
            residuals = residuals / uncert
        return residuals

    def create2dDcSurvey(self, data_type,
                         ip_type=None, pp=None,
                         no_remote=False, line_dir=None, doff=0):
        """
        Loads a dias data file to a SimPEG "srcList" class

        Input:
        datatype = Choose either IP or DC

        Note: elevation is +ve for simPEG inversion

        """
        # doff = 0                                    # in case offset is require
        srcLists = []                               # Injections + dipoles
        data = []                                   # data from file
        d_weights = []                              # weights for the data
        dc_voltage_for_ip = []
        num_rdg = len(self.readings)
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            num_dipole_count = 0
            for i in range(num_dipole):
                if data_type == "DC":
                    if (self.readings[k].Vdp[i].flagRho == "Accept"):
                        num_dipole_count += 1
                if data_type == "IP":
                    if self.readings[k].Vdp[i].flagMx == "Accept":
                        num_dipole_count += 1
            rx = np.zeros((num_dipole_count, 4))
            if line_dir.lower() == 'ns':
                tx = np.array([self.readings[k].Idp.Tx1y,
                               self.readings[k].Idp.Tx1Elev - doff,
                               self.readings[k].Idp.Tx2y,
                               self.readings[k].Idp.Tx2Elev - doff])
            elif line_dir.lower() == 'ew':
                tx = np.array([self.readings[k].Idp.Tx1x,
                               self.readings[k].Idp.Tx1Elev - doff,
                               self.readings[k].Idp.Tx2x,
                               self.readings[k].Idp.Tx2Elev - doff])
            else:
                print('[ERROR] - line directions are ew of ns')
                return
            cnt = 0
            for i in range(num_dipole):
                if data_type == "DC":
                    if (self.readings[k].Vdp[i].flagRho == "Accept"):
                        if line_dir.lower() == 'ns':
                            rx[cnt, :] = [self.readings[k].Vdp[i].Rx1y,
                                          # 0.,
                                          self.readings[k].Vdp[i].Rx1Elev - doff,
                                          self.readings[k].Vdp[i].Rx2y,
                                          # 0.,                                          
                                          self.readings[k].Vdp[i].Rx2Elev - doff]
                        elif line_dir.lower() == 'ew':
                            rx[cnt, :] = [self.readings[k].Vdp[i].Rx1x,
                                          # 0.,                            
                                          self.readings[k].Vdp[i].Rx1Elev - doff,
                                          self.readings[k].Vdp[i].Rx2x,
                                          # 0.,                                          
                                          self.readings[k].Vdp[i].Rx2Elev - doff]
                        else:
                            print('[ERROR] - line directions are ew of ns')
                            return
                        Vp = self.readings[k].Vdp[i].Vp
                        data.append(Vp / self.readings[k].Vdp[i].In)
                        d_weights.append((self.readings[k].Vdp[i].Vp_err +
                                          self.readings[k].Vdp[i].In_err) / 100.)
                        cnt += 1
                if data_type == "IP":
                    if self.readings[k].Vdp[i].flagMx == "Accept":
                        if line_dir.lower() == 'ns':
                            rx[cnt, :] = [self.readings[k].Vdp[i].Rx1y,
                                          # 0.,                            
                                          self.readings[k].Vdp[i].Rx1Elev - doff,
                                          self.readings[k].Vdp[i].Rx2y,
                                          # 0.,                                          
                                          self.readings[k].Vdp[i].Rx2Elev - doff]
                        elif line_dir.lower() == 'ew':
                            rx[cnt, :] = [self.readings[k].Vdp[i].Rx1x,
                                          # 0.,                            
                                          self.readings[k].Vdp[i].Rx1Elev - doff,
                                          self.readings[k].Vdp[i].Rx2x,
                                          # 0.,                                          
                                          self.readings[k].Vdp[i].Rx2Elev - doff]
                        else:
                            print('[ERROR] - line directions are ew of ns')
                            return
                        if ip_type is None:
                            data.append(self.readings[k].Vdp[i].Mx)
                            Vp = self.readings[k].Vdp[i].Vp * 1e-3
                            In = self.readings[k].Vdp[i].In * 1e-3
                            dc_voltage_for_ip.append(Vp / In)
                        elif ip_type == "decay":
                            In = self.readings[k].Vdp[i].In * 1e-3
                            Vp = self.readings[k].Vdp[i].Vp * 1e-3
                            data.append(self.readings[k].Vdp[i].Vs / Vp)
                            dc_voltage_for_ip.append(Vp / In)
                        elif ip_type == "voltage":
                            Vp = self.readings[k].Vdp[i].Vp * 1e-3
                            In = self.readings[k].Vdp[i].In * 1e-3
                            data.append((self.readings[k].Vdp[i].Mx * 1e-3) *
                                        (Vp / In))
                            dc_voltage_for_ip.append(Vp / In)

                        d_weights.append((self.readings[k].Vdp[i].Mx *
                                          (self.readings[k].Vdp[i].Mx_err /
                                           100.0)) / 1e3)
                        cnt += 1
            # check if remotes are required to be removed
            # print(rx.shape, tx.shape)
            if rx.shape[0] != 0:
                if pp is None:
                    Rx = DC.Rx.Dipole_ky(rx[:, :2], rx[:, 2:])    # create dipole list
                    if no_remote:
                        srcLists.append(DC.Src.Pole([Rx], tx[:2]))
                    else:
                        srcLists.append(DC.Src.Dipole([Rx], tx[:2], tx[2:]))
                else:
                    cnt_pp = 0
                    cnt_dp = 0
                    for idjj in range(rx[:, 0].size):
                        if np.all(rx[idjj, 2:] == pp):
                            cnt_pp += 1
                        else:
                            cnt_dp += 1

                    rx_pole = np.zeros((cnt_pp, 4))
                    rx_dp = np.zeros((cnt_dp, 4))
                    cnt_pp = 0
                    cnt_dp = 0
                    for idj in range(rx[:, 0].size):
                        if np.all(rx[idj, 2:] == pp):
                            rx_pole[cnt_pp, :] = rx[idj, :]
                            cnt_pp += 1
                        else:
                            rx_dp[cnt_dp, :] = rx[idj, :]
                            cnt_dp += 1
                    Rx_pp = DC.Rx.Pole_ky(rx_pole[:, :2])
                    Rx_dp = DC.Rx.Dipole_ky(rx_dp[:, :2], rx_dp[:, 2:])
                   
                    if no_remote:
                        srcLists.append(DC.Src.Pole([Rx_dp], tx[:2]))
                        srcLists.append(DC.Src.Pole([Rx_pp], tx[:2]))
                    else:
                        srcLists.append(DC.Src.Dipole([Rx_dp], tx[:2], tx[2:]))
                        srcLists.append(DC.Src.Dipole([Rx_pp], tx[:2], tx[2:]))
            else:
                print('[INFO] found bad IP reading')

        survey = DC.Survey_ky(srcLists)          # creates the survey
        # check if data is IP
        if data_type == "IP":
            survey = IP.from_dc_to_ip_survey(survey, dim="2.5D")
        survey.dobs = np.float64(np.asarray(data))                 # assigns data
        if ip_type == "decay":
            survey.dc_voltage = np.float64(np.asarray(dc_voltage_for_ip))
        else:
            survey.dc_voltage = np.float64(np.asarray(dc_voltage_for_ip))
        survey.std = np.asarray(d_weights)             # assign data weights
        survey.eps = 0.001

        return survey

    def createDcSurvey(self, data_type, ip_type=None, pp=None, no_remote=False, reject=None):
        """
        Loads a dias data file to a SimPEG "srcList" class

        Input:
        datatype = Choose either IP or DC

        Note: elevation is +ve for simPEG inversion

        """
        doff = 0                                    # in case offset is require
        srcLists = []                               # Injections + dipoles
        data = []                                   # data from file
        dc_voltage_for_ip = []
        d_weights = []                              # weights for the data
        num_rdg = len(self.readings)
        minE = self.readings[0].Idp.Tx2East
        minN = self.readings[0].Idp.Tx2North
        maxN = minN
        maxE = minE
        if reject is not None:
            print('[INFO] Dc is using Mx as flag')
        for k in range(num_rdg):
            num_dipole = len(self.readings[k].Vdp)
            num_dipole_count = 0
            for i in range(num_dipole):
                if data_type == "DC":
                    if reject is None:
                        if (self.readings[k].Vdp[i].flagRho == "Accept"):
                            num_dipole_count += 1
                    elif reject == "app_mx":
                        if self.readings[k].Vdp[i].flagMx == "Accept":
                            num_dipole_count += 1
                if data_type == "IP":
                    if self.readings[k].Vdp[i].flagMx == "Accept":
                        num_dipole_count += 1
            rx = np.zeros((num_dipole_count, 6))
            tx = np.array([self.readings[k].Idp.Tx1East,
                           self.readings[k].Idp.Tx1North,
                           self.readings[k].Idp.Tx1Elev - doff,
                           self.readings[k].Idp.Tx2East,
                           self.readings[k].Idp.Tx2North,
                           self.readings[k].Idp.Tx2Elev - doff])
            if self.readings[k].Idp.Tx1East > maxE:
                maxE = self.readings[k].Idp.Tx1East
            if self.readings[k].Idp.Tx1East < minE:
                minE = self.readings[k].Idp.Tx1East
            if self.readings[k].Idp.Tx1North > maxN:
                maxN = self.readings[k].Idp.Tx1North
            if self.readings[k].Idp.Tx1North < minN:
                minN = self.readings[k].Idp.Tx1North
            cnt = 0
            for i in range(num_dipole):
                if data_type == "DC":
                    if reject is None:
                        if (self.readings[k].Vdp[i].flagRho == "Accept"):
                            rx[cnt, :] = [self.readings[k].Vdp[i].Rx1East,
                                          self.readings[k].Vdp[i].Rx1North,
                                          self.readings[k].Vdp[i].Rx1Elev - doff,
                                          self.readings[k].Vdp[i].Rx2East,
                                          self.readings[k].Vdp[i].Rx2North,
                                          self.readings[k].Vdp[i].Rx2Elev - doff]
                            Vp = self.readings[k].Vdp[i].Vp
                            data.append(Vp / self.readings[k].Vdp[i].In)
                            d_weights.append((self.readings[k].Vdp[i].Vp_err +
                                              self.readings[k].Vdp[i].In_err) / 100.)
                            cnt += 1
                    elif reject == "app_mx":
                        if (self.readings[k].Vdp[i].flagMx == "Accept"):
                            rx[cnt, :] = [self.readings[k].Vdp[i].Rx1East,
                                          self.readings[k].Vdp[i].Rx1North,
                                          self.readings[k].Vdp[i].Rx1Elev - doff,
                                          self.readings[k].Vdp[i].Rx2East,
                                          self.readings[k].Vdp[i].Rx2North,
                                          self.readings[k].Vdp[i].Rx2Elev - doff]
                            Vp = self.readings[k].Vdp[i].Vp
                            data.append(Vp / self.readings[k].Vdp[i].In)
                            d_weights.append((self.readings[k].Vdp[i].Vp_err +
                                              self.readings[k].Vdp[i].In_err) / 100.)
                            cnt += 1
                if data_type == "IP":
                    if self.readings[k].Vdp[i].flagMx == "Accept":
                        rx[cnt, :] = [self.readings[k].Vdp[i].Rx1East,
                                      self.readings[k].Vdp[i].Rx1North,
                                      self.readings[k].Vdp[i].Rx1Elev - doff,
                                      self.readings[k].Vdp[i].Rx2East,
                                      self.readings[k].Vdp[i].Rx2North,
                                      self.readings[k].Vdp[i].Rx2Elev - doff]
                        if ip_type is None:
                            data.append(self.readings[k].Vdp[i].Mx)
                            Vp = self.readings[k].Vdp[i].Vp * 1e-3
                            In = self.readings[k].Vdp[i].In * 1e-3
                            dc_voltage_for_ip.append(Vp / In)
                        elif ip_type == "decay":
                            Vp = self.readings[k].Vdp[i].Vp * 1e-3
                            In = self.readings[k].Vdp[i].In * 1e-3
                            data.append(self.readings[k].Vdp[i].Vs / Vp)
                            dc_voltage_for_ip.append(Vp / In)
                        elif ip_type == "voltage":
                            Vp = self.readings[k].Vdp[i].Vp * 1e-3
                            In = self.readings[k].Vdp[i].In * 1e-3
                            data.append((self.readings[k].Vdp[i].Mx * 1e-3) *
                                        (Vp / In))

                        d_weights.append((self.readings[k].Vdp[i].Mx *
                                          (self.readings[k].Vdp[i].Mx_err /
                                           100.0)) / 1e3)
                        cnt += 1
            # check if remotes are required to be removed
            if pp is None:
                Rx = DC.receivers.Dipole(rx[:, :3], rx[:, 3:])    # create dipole list
                if no_remote:
                    srcLists.append(DC.sources.Pole([Rx], tx[:3]))
                else:
                    srcLists.append(DC.sources.Dipole([Rx], tx[:3], tx[3:]))
            else:
                cnt_pp = 0
                cnt_dp = 0
                for idjj in range(rx[:, 0].size):
                    if np.all(rx[idjj, 3:] == pp):
                        cnt_pp += 1
                    else:
                        cnt_dp += 1

                rx_pole = np.zeros((cnt_pp, 6))
                rx_dp = np.zeros((cnt_dp, 6))
                cnt_pp = 0
                cnt_dp = 0
                for idj in range(rx[:, 0].size):
                    if np.all(rx[idj, 3:] == pp):
                        rx_pole[cnt_pp, :] = rx[idj, :]
                        cnt_pp += 1
                    else:
                        rx_dp[cnt_dp, :] = rx[idj, :]
                        cnt_dp += 1
                Rx_pp = DC.receivers.Pole(rx_pole[:, :3])
                Rx_dp = DC.receivers.Dipole(rx_dp[:, :3], rx_dp[:, 3:])
                
                if no_remote:
                    srcLists.append(DC.sources.Pole([Rx_dp], tx[:3]))
                    srcLists.append(DC.sources.Pole([Rx_pp], tx[:3]))
                else:
                    srcLists.append(DC.sources.Dipole([Rx_dp], tx[:3], tx[3:]))
                    srcLists.append(DC.sources.Dipole([Rx_pp], tx[:3], tx[3:]))

        survey = DC.Survey(srcLists)          # creates the survey
        # check if data is IP
        if data_type == "IP":
            survey = IP.from_dc_to_ip_survey(survey, dim="3D")
        survey.dobs = np.float64(np.asarray(data))                 # assigns data
        if ip_type == "decay":
            survey.dc_voltage = np.float64(np.asarray(dc_voltage_for_ip))
        if ip_type == None:
            survey.dc_voltage = np.float64(np.asarray(dc_voltage_for_ip))
        survey.std = np.asarray(d_weights)             # assign data weights
        survey.eps = 0.001

        return survey

    def plotGpsOverDatabaseLocations(self, gps_input=None, rx=True, tx=True, dipole_dipole=False, local=False, line_dir=None):
        if gps_input is not None:
            if rx and tx:
                tx = self.getSources(local=local)
                rx = self.getDipoles(local=local)
                plt.plot(tx[:, 0], tx[:, 1], 'oc')
                if dipole_dipole:
                    plt.plot(tx[:, 3], tx[:, 4], 'oc')
                plt.plot(rx[:, 0], rx[:, 1], 'oc')
                plt.plot(gps_input[0, :], gps_input[1, :], '+k')
                plt.xlabel('Easting (m)')
                plt.ylabel('Northing (m)')
                plt.title('Database Locations & GPS View')
                plt.grid()
                plt.show()
            elif rx:
                rx = self.getDipoles(local=local)
                plt.plot(rx[:, 0], rx[:, 1], 'oc')
                plt.plot(gps_input[0, :], gps_input[1, :], '+k')
                plt.xlabel('Easting (m)')
                plt.ylabel('Northing (m)')
                plt.title('Database Rx Locations & GPS View')
                plt.grid()
                plt.show()
            elif tx:
                tx = self.getSources(local=local)
                plt.plot(tx[:, 0], tx[:, 1], 'om')
                if dipole_dipole:
                    plt.plot(tx[:, 3], tx[:, 4], 'om')
                plt.plot(gps_input[0, :], gps_input[1, :], '+k')
                plt.xlabel('Easting (m)')
                plt.ylabel('Northing (m)')
                plt.grid()
                plt.title('Database Tx Locations & GPS View')
                plt.show()
        else:
            print("[INFO] Please provide gps input data!!!!")

    def compareInjectedCurrents2fieldDDN(self, ddn=None, reject=None, line_dir='ew'):
        if ddn is None:
            print("[INFO] Please provide a DDN filed from the field")
        else:
            # make a dictionary of the DDN for easy look up
            record_dict = {}
            for row in range(ddn.shape[0]):
                record_dict[str(int(ddn[row, 0]))] = ddn[row, 3]

            # create list of difference data
            diff_in = []
            in_ = []
            in_local = []
            rec_local = []
            diff_recs = []
            num_rdg = len(self.readings)
            for rdg in range(num_rdg):
                num_dipole = len(self.readings[rdg].Vdp)

                if line_dir == 'ew':
                    in_local.append(self.readings[rdg].Idp.Tx1x)
                    rec_local.append(int(self.readings[rdg].MemNumber))
                else:
                    in_local.append(self.readings[rdg].Idp.Tx1y)
                    rec_local.append(int(self.readings[rdg].MemNumber))

                for dp in range(num_dipole):
                    if reject == 'app_rho':
                        if self.readings[rdg].Vdp[dp].flagRho == 'Accept':
                            in_database = self.readings[rdg].Vdp[dp].In
                            dict_id = str(self.readings[rdg].MemNumber)
                            in_ddn = record_dict[dict_id]
                            diff_in.append(in_database - in_ddn)
                            diff_recs.append(self.readings[rdg].MemNumber)
                            in_.append(in_database)
                    if reject == 'app_mx':
                        if self.readings[rdg].Vdp[dp].flagMx == 'Accept':
                            in_database = self.readings[rdg].Vdp[dp].In
                            dict_id = str(self.readings[rdg].MemNumber)
                            in_ddn = record_dict[dict_id]
                            diff_in.append(in_database - in_ddn)
                            in_.append(in_database)
                            diff_recs.append(int(self.readings[rdg].MemNumber))
                    else:
                        in_database = np.abs(self.readings[rdg].Vdp[dp].In)
                        dict_id = str(self.readings[rdg].MemNumber)
                        try:
                            in_ddn = record_dict[dict_id]
                            diff_in.append(np.abs(in_database - in_ddn) / ((in_database + in_ddn) / 2) * 100)
                            in_.append(in_database)
                            diff_recs.append(int(self.readings[rdg].MemNumber))
                        except KeyError:
                            print('No matching mem for: {0}'.format(dict_id))

            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            ax1.plot(diff_recs, diff_in, '*g')
            ax1.set_title(r'$\Delta$' + 'I || mean: ' + str(np.mean(diff_in)))
            ax1.grid()
            ax1.set_ylabel('percent difference')
            # ax1.set_ylim([np.min(in_), np.max(in_)])
            ax2.plot(diff_recs, in_, '*')
            ax2.plot(ddn[:, 0], ddn[:, 3], 'om')
            ax2.set_title('I_ddn & I_db')
            ax2.legend(['db', 'ddn'])
            ax2.set_xlabel('record #')
            ax2.set_ylabel('In (mA)')
            ax2.grid()
            plt.show()

            plt.plot(rec_local, in_local, '^-')
            plt.plot(ddn[:, 0], ddn[:, 2], '+-r')
            plt.xlabel('Record #')
            plt.ylabel('Station')
            plt.grid()
            plt.title('Record NUmber Vs Station')
            plt.show()

    def plotLabeledInjections(self, gps=None):
        if gps is not None:
            for rdg in range(len(self.readings)):
                tx1x = self.readings[rdg].Idp.Tx1East
                tx1y = self.readings[rdg].Idp.Tx1North
                rec_id = str(int(self.readings[rdg].MemNumber))
                stn_x_id = str(int(self.readings[rdg].Idp.Tx1x))
                stn_y_id = str(int(self.readings[rdg].Idp.Tx1y))
                id_ = rec_id + "-" + stn_x_id + "-" + stn_y_id
                plt.plot(tx1x, tx1y, '-dg')
                plt.text(tx1x, tx1y, '%s' % id_)
                plt.axis('equal')
                plt.title("Injection ID map")
                plt.xlabel('Easting (m)')
                plt.ylabel('Northing (m)')
                # plt.rcParams.update({'font.size': 4})
            plt.show()
        else:
            print("to do")

    def checkDipoleAppRhoPolarityPerReading(self, num_readings=None,
                                            gps_locations=None,
                                            dipole_dipole=False):
        if num_readings is not None:
            num_rdg = num_readings
            for rdg in range(int(num_rdg)):
                self.readings[rdg].createNodeDB()
                plt.figure(figsize=(9, 8))
                plt.axis('equal')
                for idx in range(len(self.readings[rdg].node_db)):
                    x = self.readings[rdg].node_locs[idx][0]
                    y = self.readings[rdg].node_locs[idx][1]
                    plt.plot(x, y, 'k*')
                    plt.text(x, y, '%s' % self.readings[rdg].node_db[idx])
                plt.title("rec: {0}".format(self.readings[rdg].MemNumber))
                for dp in range(len(self.readings[rdg].Vdp)):
                    rx1x = self.readings[rdg].Vdp[dp].Rx1East
                    rx1y = self.readings[rdg].Vdp[dp].Rx1North
                    rx2x = self.readings[rdg].Vdp[dp].Rx2East
                    rx2y = self.readings[rdg].Vdp[dp].Rx2North
                    In = self.readings[rdg].Idp
                    if self.readings[rdg].Vdp[dp].calcRho(In) > 0:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-ob')
                    else:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-or')
                if gps_locations is not None:
                    plt.plot(gps_locations[0, :], gps_locations[1, :], '+k')
                if dipole_dipole:
                    tx1x = self.readings[rdg].Idp.Tx1East
                    tx1y = self.readings[rdg].Idp.Tx1North
                    tx2x = self.readings[rdg].Idp.Tx2East
                    tx2y = self.readings[rdg].Idp.Tx2North
                    plt.plot([tx1x, tx2x], [tx1y, tx2y], '-dk')
                else:
                    tx1x = self.readings[rdg].Idp.Tx1East
                    tx1y = self.readings[rdg].Idp.Tx1North
                    plt.plot(tx1x, tx1y, 'dk')
                plt.show()
        else:
            num_rdg = len(self.readings)
            for rdg in range(int(num_rdg)):
                self.readings[rdg].createNodeDB()
                for idx in range(len(self.readings[rdg].node_db)):
                    x = self.readings[rdg].node_locs[idx][0]
                    y = self.readings[rdg].node_locs[idx][1]
                    plt.plot(x, y, 'k*')
                    plt.text(x, y, '%s' % self.readings[rdg].node_db[idx])
                plt.title("rec: {0}".format(self.readings[rdg].MemNumber))
                for dp in range(len(self.readings[rdg].Vdp)):
                    rx1x = self.readings[rdg].Vdp[dp].Rx1East
                    rx1y = self.readings[rdg].Vdp[dp].Rx1North
                    rx2x = self.readings[rdg].Vdp[dp].Rx2East
                    rx2y = self.readings[rdg].Vdp[dp].Rx2North
                    In = self.readings[rdg].Idp
                    if self.readings[rdg].Vdp[dp].calcRho(In) > 0:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-ob')
                    else:
                        plt.plot([rx1x, rx2x], [rx1y, rx2y], '-or')
                if gps_locations is not None:
                    plt.plot(gps_locations[0, :], gps_locations[1, :], '+k')
                if dipole_dipole:
                    tx1x = self.readings[rdg].Idp.Tx1East
                    tx1y = self.readings[rdg].Idp.Tx1North
                    tx2x = self.readings[rdg].Idp.Tx2East
                    tx2y = self.readings[rdg].Idp.Tx2North
                    plt.plot([tx1x, tx2x], [tx1y, tx2y], '-dk')
                else:
                    tx1x = self.readings[rdg].Idp.Tx1East
                    tx1y = self.readings[rdg].Idp.Tx1North
                    plt.plot(tx1x, tx1y, 'dk')
                plt.show()

    def plotHistogramOfDataVitals(self, reject='app_rho',
                                  log_rho=False, log_vp=False):
        fig = plt.figure()
        ax1 = fig.add_subplot(241)
        ax2 = fig.add_subplot(242)
        ax3 = fig.add_subplot(243)
        ax4 = fig.add_subplot(244)
        ax5 = fig.add_subplot(245)
        ax6 = fig.add_subplot(246)
        ax7 = fig.add_subplot(247)
        ax8 = fig.add_subplot(248)

        ax1.set_title(r'$\rho$')
        ax2.set_title("error V/I")
        ax3.set_title("K")
        ax4.set_title("V")
        ax5.set_title(r'$\eta$')
        ax6.set_title("Null")
        ax7.set_title("I")
        ax8.set_title("a")

        rho = self.getApparentResistivity(reject=reject)
        vp_i_error = self.getVpDivInErrors(reject=reject)
        gk = self.getGeometricFactor(reject=reject)
        voltages = self.getVoltages(reject=reject)
        mx = self.getApparentChageability(reject=reject)
        null = self.getNullFactors(reject=reject)
        current = self.getCurrents(reject=reject)
        a = self.getAspacings(reject=reject)

        if log_rho:
            ax1.hist(np.log(np.abs(rho)), 50)
            ax1.set_title("log(" + r'$\rho$' + ")")
        else:
            ax1.hist(rho, 50)
        ax2.hist(vp_i_error, 25)
        ax3.hist(gk, 50)
        if log_vp:
            ax4.hist(np.log(voltages), 50)
            ax4.set_title("log(V)")
        else:
            ax4.hist(voltages, 50)
        ax5.hist(mx, range=None)
        ax6.hist(null, 50)
        ax7.hist(current, 50)
        ax8.hist(a)

        plt.show()

    def compareRecords(self, rec1, rec2,
                       histogram=True,
                       histogram_ratio=True,
                       histo_bins_dc=30,
                       histo_bins_ip=30,
                       histo_bins_ratio=30,
                       histo_bins_xc=30,
                       ylimit_dc=None,
                       ylimit_ip=None):
        assert type(rec1) == str, "input must be string"
        # lets find the index of these reading
        id1 = -1
        id2 = -1
        # find the indexes of requested readings
        for idx in range(len(self.readings)):
            print(self.readings[idx].MemNumber, rec1)
            if self.readings[idx].MemNumber == rec1:
                id1 = idx
            if self.readings[idx].MemNumber == rec2:
                id2 = idx
        print(id1, id2)
        # check that the records were found
        if id1 != -1 and id2 != -1:
            # set plot
            fig = plt.figure()
            fig.suptitle('Readings: ' + rec1 + ' (blue) & ' + rec2 +
                         ' (red)', fontsize=16)
            ax1 = fig.add_subplot(231)
            ax2 = fig.add_subplot(232)
            ax3 = fig.add_subplot(233)
            ax4 = fig.add_subplot(234)
            ax5 = fig.add_subplot(235)
            ax6 = fig.add_subplot(236)
            # plot all decays from one
            for idx in range(len(self.readings[id1].Vdp)):
                ax1.plot(self.window_center, self.readings[id1].Vdp[idx].Vs /
                         self.readings[id1].Vdp[idx].Vp, '-ob')
            for idx in range(len(self.readings[id2].Vdp)):
                ax1.plot(self.window_center, self.readings[id2].Vdp[idx].Vs /
                         self.readings[id2].Vdp[idx].Vp, '-or')
            ax1.set_title('Decay')
            ax1.set_xlabel('time (ms)')
            ax1.set_ylabel('Voltage (mV/V)')

            # now Rho compare
            dc1 = []
            dc2 = []
            for idx in range(len(self.readings[id1].Vdp)):
                dc1.append(self.readings[id1].Vdp[idx].calcRho(self.readings[id1].Idp))
            for idx in range(len(self.readings[id2].Vdp)):
                dc2.append(self.readings[id2].Vdp[idx].calcRho(self.readings[id2].Idp))
            if histogram:
                ax2.hist(dc1, histo_bins_dc)
                ax2.hist(dc2, histo_bins_dc,  fc=[1, 0, 0, 0.5])
            else:
                ax2.plot(dc1, '*b')
                ax2.plot(dc2, 'or')

            ax2.set_title('Resistivity compare')
            ax2.set_ylabel('ohm-m')
            if ylimit_dc is not None:
                ax2.set_ylim(ylimit_dc)

            # now Mx
            ip1 = []
            ip2 = []
            for idx in range(len(self.readings[id1].Vdp)):
                ip1.append(self.readings[id1].Vdp[idx].Mx)
            for idx in range(len(self.readings[id2].Vdp)):
                ip2.append(self.readings[id2].Vdp[idx].Mx)
            
            if histogram:
                ax3.hist(ip1, histo_bins_ip)
                ax3.hist(ip2, histo_bins_ip,  fc=[1, 0, 0, 0.5])
            else:
                ax3.plot(ip1, '*b')
                ax3.plot(ip2, 'or')
            ax3.set_title('Mx compare')
            ax3.set_ylabel('mV/V')
            if ylimit_ip is not None:
                ax3.set_ylim(ylimit_ip)

            # now subtract DC
            subtracted_values = []
            for idx in range(len(self.readings[id1].Vdp)):
                node_id1 = self.readings[id1].Vdp[idx].Rx1File[:2]
                node_id2 = self.readings[id1].Vdp[idx].Rx2File[:2]
                for idj in range(len(self.readings[id2].Vdp)):
                    node_id1_ = self.readings[id2].Vdp[idj].Rx1File[:2]
                    node_id2_ = self.readings[id2].Vdp[idj].Rx2File[:2]
                    if node_id1_ == node_id1 and node_id2_ == node_id2:
                        # print('subtracting {0} - {1} & {2} - {3}'.format(node_id1, node_id1_, node_id2, node_id2_))
                        diff_value = (self.readings[id1].Vdp[idx].calcRho(self.readings[id1].Idp) -
                                      self.readings[id2].Vdp[idj].calcRho(self.readings[id2].Idp))
                        value1 = self.readings[id1].Vdp[idx].calcRho(self.readings[id1].Idp)
                        value2 = self.readings[id2].Vdp[idj].calcRho(self.readings[id2].Idp)
                        try:
                            if np.abs(value1) < np.abs(value2):
                                diff_value = value1 / value2
                            else:
                                diff_value = value2 / value1
                        except ZeroDivisionError:
                            print('zero division - do you think you are not human!')
                            diff_value = 0
                        subtracted_values.append(np.abs(diff_value))

            if histogram_ratio:
                ax4.hist(subtracted_values, histo_bins_ratio)
            else:
                ax4.plot(subtracted_values, 'o')
            ax4.set_title("DC-ratio - mean:" + str(np.mean(subtracted_values)) + " median: " +
                          str(np.median(subtracted_values)))

            # now subtract IP
            subtracted_values = []
            for idx in range(len(self.readings[id1].Vdp)):
                node_id1 = self.readings[id1].Vdp[idx].Rx1File[:2]
                node_id2 = self.readings[id1].Vdp[idx].Rx2File[:2]
                for idj in range(len(self.readings[id2].Vdp)):
                    node_id1_ = self.readings[id2].Vdp[idj].Rx1File[:2]
                    node_id2_ = self.readings[id2].Vdp[idj].Rx2File[:2]
                    if node_id1_ == node_id1 and node_id2_ == node_id2:
                        # print('subtracting {0} - {1} & {2} - {3}'.format(node_id1, node_id1_, node_id2, node_id2_))
                        diff_value = (self.readings[id1].Vdp[idx].Mx -
                                      self.readings[id2].Vdp[idj].Mx)
                        value1 = np.abs(self.readings[id1].Vdp[idx].Mx)
                        value2 = np.abs(self.readings[id2].Vdp[idj].Mx)
                        try:
                            if value1 < value2:
                                diff_value = value1 / value2
                            else:
                                diff_value = value2 / value1
                            # if diff_value > 1:
                            #     print(value1, value2)
                        except ZeroDivisionError:
                            diff_value = 0
                        subtracted_values.append(np.abs(diff_value))
            
            if histogram_ratio:
                ax5.hist(subtracted_values, histo_bins_ratio)
            else:
                ax5.plot(subtracted_values, 'o')
            ax5.set_title("IP-ratio - mean:" + str(np.mean(subtracted_values)) + " median: " +
                          str(np.median(subtracted_values)))

            # now for Xcorr
            # now Rho compare
            xcor1 = []
            xcor2 = []
            for idx in range(len(self.readings[id1].Vdp)):
                xcor1.append(self.readings[id1].Vdp[idx].xcorr)
            for idx in range(len(self.readings[id2].Vdp)):
                xcor2.append(self.readings[id2].Vdp[idx].xcorr)
            ax6.hist(xcor1, histo_bins_xc)
            ax6.hist(xcor2, histo_bins_xc, fc=[1, 0, 0])
            ax6.set_title('Cross-Correlation')
            plt.show()
        else:
            print('[WARNING] COuld not find requested record IDs')


class Jreadtxtline:
    """
    Class specifically for reading a line of text from a dias file

    """

    def __init__(self, hdrLine, dataLine):
        # make a structure with the header inputs
        self.Vs = []
        for n in range(len(hdrLine)):
            if hdrLine[n].find("Vs") == 0:
                try:
                    self.Vs.append(float(dataLine[n]))
                except ValueError:
                    self.Vs.append(float(-99.9))

            else:
                setattr(self, hdrLine[n], dataLine[n])
