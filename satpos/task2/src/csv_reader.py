import csv
import numpy as np
import matplotlib.pyplot as plt


def csvl(fpath):
    years = list()
    val = list()
    error = list()

    with open(fpath, "r") as f:
        samples = 0
        samples_per_day = 86400./300.
        for l in csv.reader(f, skipinitialspace=True, delimiter=" "):
            samples += 1
            years.append(float(l[0]))
            val.append(float(l[2]))
            error.append(float(l[3]))
    return [np.array(years, dtype=np.float64),
            np.array(val, dtype=np.float64),
            np.array(error, dtype=np.float64),
            samples/samples_per_day]


class station(object):
    def __init__(self, h_correction, filepath):
        [years, val, error, days] = csvl(filepath)
        self.years = years
        self.val = val + h_correction
        self.error = error
        self.n_days = days

    def fit_data(self, days, optimize):
        end = self.idx_from_day(days)
        t = self.years[:end]
        v = self.val[:end]
        e = self.error[:end] * np.sqrt(1+(end-np.arange(end)))
        c_bf = np.polyfit(t, v, 0, w=1.0/e, full=False)
        if optimize:
            y_bf = self.__data_from_coeff(t, c_bf)
            for i in range(0, 20):
                c = np.polynomial.polynomial.polyfit(t, v, i, w=1.0/e, full=False)
                y = self.__data_from_coeff(t, c)
                if (np.linalg.norm(y - v) < np.linalg.norm(y_bf - v)):
                    y_bf = y
                    c_bf = c
        self.valfit = self.__data_from_coeff(self.years, c_bf)
        return self.valfit

    def __data_from_coeff(self, t, c):
        y = np.zeros(np.size(t, axis=0))
        t_pow = np.ones(np.size(t, axis=0))
        for i in c:
            y += i*t_pow
            t_pow *= t
        return y

    def idx_from_day(self, day):
        return int(day/self.n_days*np.size(self.years, axis=0))


def do_station(stn_obj, ax):
    last_day = stn_obj.n_days-1
    i = stn_obj.idx_from_day(last_day)
    stn_obj.fit_data(last_day, True)
    ax.errorbar(stn_obj.years[:i], stn_obj.val[:i], yerr=stn_obj.error[:i],
                label="GPS data day 1-6",
                linestyle="-", ecolor="#0000ff")
    ax.plot(stn_obj.years, stn_obj.valfit, color="#ff0000", linestyle="--",
            label="Polyfit %.0f/%.0f days" % (last_day, stn_obj.n_days))
    stn_obj.fit_data(stn_obj.n_days, True)
    ax.plot(stn_obj.years, stn_obj.valfit, color="#00ff00", linestyle="-.",
            label="Polyfit %.0f/%.0f days" % (stn_obj.n_days, stn_obj.n_days))
    ax.errorbar(stn_obj.years[i:], stn_obj.val[i:], yerr=stn_obj.error[i:],
                label="GPS data day 7",
                linestyle="-", ecolor="#000000", color="#000000")
    ax.legend()

    
if __name__ == "__main__":
    fig = plt.figure()
    axis_reso = fig.add_subplot(211)
    axis_reso.set_title("RESO")
    axis_mdo1 = fig.add_subplot(212)
    axis_mdo1.set_title("MDO1")
    reso = station(2.2941, "../res/trop/RESO_trop")
    do_station(reso, axis_reso)
    mdo1 = station(1.8223, "../res/trop/MDO1_trop")
    do_station(mdo1, axis_mdo1)
    plt.show()
