from imports import *

def import_all_flares(fname):
    df = pd.read_csv(fname)
    return df

def create_flaredict(amp, fwhm, istart, istop, energy, modenergy, sat):
    flaredict = {"amp": amp, "fwhm": fwhm, 'istart': istart, 'istop': istop, 'energy': energy,
                 'model_energy': modenergy, 'sat': sat}
    return flaredict

def write_to_file(tpeaks, amps, fwhms, energies, wn, rn, trec, arec, frec, erec, recovered, svname):
    if os.path.exists(svname):
        print(svname, "exists!")
        with open(svname, 'a') as fd:
            writer = csv.writer(fd)
            for i in range(len(tpeaks)):
                csvrow = [tpeaks[i], fwhms[i], amps[i], energies[i], recovered[i], trec[i], arec[i], frec[i], erec[i],
                          wn[i], rn[i]]
                writer.writerow(csvrow)
        fd.close()

    else:
        print(svname, " doesn't exist!")
        with open(svname, 'w') as fd:
            writer = csv.writer(fd)
            writer.writerow(
                ['TPEAK-2458000', 'FWHM_INJ', 'AMP_INJ', "ENERGY_INJ", "RECOVERED", 'TPEAK_REC', 'FWHM_REC', 'AMP_REC',
                 'ENERGY_REC', 'WN', 'RN'])
            for i in range(len(tpeaks)):
                csvrow = [tpeaks[i], fwhms[i], amps[i], energies[i], recovered[i], trec[i], arec[i], frec[i], erec[i],
                          wn[i], rn[i]]
                print(csvrow)
                writer.writerow(csvrow)
        fd.close()


def write_to_file_recovered(tpeaks, amps, fwhms, energies, svname):
    if os.path.exists(svname):
        print(svname, "exists!")
        with open(svname, 'a') as fd:
            writer = csv.writer(fd)
            for i in range(len(tpeaks)):
                csvrow = [tpeaks[i], fwhms[i], amps[i], energies[i]]
                writer.writerow(csvrow)
        fd.close()

    else:
        print(svname, " doesn't exist!")
        with open(svname, 'w') as fd:
            writer = csv.writer(fd)
            writer.writerow(['TPEAK-2458000', 'FWHM', 'AMP', "ENERGY"])
            for i in range(len(tpeaks)):
                csvrow = [tpeaks[i], fwhms[i], amps[i], energies[i]]
                print(csvrow)
                writer.writerow(csvrow)
        fd.close()