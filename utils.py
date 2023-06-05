from imports import *


def bin_data(t, y, b):
    """
    Function to bin a variable, y, based on binning time series, t, with bin size, b.
    :param t: time series
    :param y: flux series
    :param b: bin size in minutes
    :return: binned time, binned flux and error in bins
    """
    mins_jd = float(b) / 1440.
    t = np.array(t)
    y = np.array(y)

    split = []
    sorted_t = t
    sorted_y = y
    start = sorted_t[0]
    nextbin = sorted_t[np.absolute(sorted_t - start) > mins_jd]

    while nextbin != []:
        start = start + mins_jd
        ind_st = np.argmax(sorted_t > start)
        if len(split) > 0:
            if ind_st != split[-1]:
                split.append(ind_st)
                time = sorted_t[ind_st:]
        else:
            split.append(ind_st)
            time = sorted_t[ind_st:]
        nextbin = time[np.absolute(time - start) > mins_jd]

    times = np.split(sorted_t, split)
    ys = np.split(sorted_y, split)

    bins = np.zeros(len(times))
    binned_y = np.zeros(len(times))
    binned_err = np.zeros(len(times))

    for i in range(len(times)):
        if len(ys[i]) > 0:
            try:
                bins[i] = np.nanmedian(times[i])
                binned_y[i] = np.nanmedian(ys[i])
                n = len(times[i])
                # error in median
                binned_err[i] = 1.253 * np.nanstd(ys[i]) / np.sqrt(n)
            except Exception as e:
                pass

    bin_t = bins[binned_y != 0]
    bin_e = binned_err[binned_y != 0]
    bin_y = binned_y[binned_y != 0]

    return bin_t, bin_y, bin_e

def in_bins(bins, x, y, width):
    bindict = {}

    for b in bins:
        bindict[b] = []
        for i in range(len(y)):
            if x[i] > b - (0.5 * width) and x[i] < b + (0.5 * width):
                if b in bindict.keys():
                    bindict[b].append(y[i])
                else:
                    bindict[b] = [y[i]]

    return bindict


def split_multilc(t):
    start = t[0]
    t = np.array(t)
    nextdays = t[np.absolute(t - start) > 0.5]
    split = []
    while nextdays != []:
        start = nextdays[0]
        ind_st = np.where(t == start)[0][0]
        split.append(ind_st)
        time = t[ind_st:]
        nextdays = time[np.absolute(time - start) > 0.5]

    times = np.split(t, split)
    # xs=np.split(x,split)

    return times, split

def get_total_observed_time(t):
    ts, split = split_multilc(t)
    total_time = 0
    for day in ts:
        t_range = day[-1] - day[0]
        total_time = total_time + t_range

    return total_time

def find_nearest(array, value):
    """
    Find the closest value in an array to a given value.
    :param array: array to search for the nearest value
    :param value: the value to find the closest to in array
    :return: the array index of the item closest to the passed value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def clipped_std(x):
    return np.ma.std(sigma_clip(x, sigma=3))

def running_box(x, y, boxsize, operation):
    if operation == "std":
        op = np.nanstd
    elif operation == 'median':
        op = np.nanmedian
    elif operation == 'mean':
        op = np.nanmean
    elif operation == "clipped_std":
        op = clipped_std
    else:
        print("No running function selected - choose std, median or mean.")
        return np.nan

    l = len(x)
    boxsize = boxsize / 2
    dy = np.zeros((np.size(x)))

    for i in range(0, l):
        s = x[i]
        box_s = find_nearest(x, s - boxsize)
        box_e = find_nearest(x, s + boxsize)

        # calculate [OPERATION] for the current window
        try:
            y_ext = y[box_s:box_e + 1]
            d = op(y_ext)
            if not np.isnan(d):
                dy[i] = d

        except Exception as e:
            print(e)

    dy = np.ma.masked_where(dy == 0, dy)
    # dy = np.ma.masked_where(np.isnan(dy), dy)
    dy = np.ma.masked_invalid(dy)
    # dy = dy.filled(np.nanmedian(dy))
    # print("RUNNING " + operation)
    # print("FILL VALUE = " + str(np.nanmedian(dy)))
    return dy


def calculate_running_rms(x, y, boxsize):
    dy = running_box(x, y, boxsize, 'std')
    return dy


def calculate_running_clipped_rms(x, y, boxsize):
    dy = running_box(x, y, boxsize, 'clipped_std')
    return dy


def calculate_running_median(x, y, boxsize):
    dy = running_box(x, y, boxsize, 'median')
    return dy


def calculate_running_mean(x, y, boxsize):
    dy = running_box(x, y, boxsize, 'mean')
    return dy

def get_response(fname):
    lam, response = [], []
    with open(fname, 'r') as rfile:
        reader = csv.reader(rfile)
        for row in reader:
            try:
                # print(float(row[0]))
                lam.append(0.000001 * float(row[0]))  # m
                response.append(float(row[1]))
            except:
                pass
    rfile.close()
    return np.array(lam), np.array(response)

def tel_response(lam, response, wavelength):
    ind = find_nearest(lam, wavelength)
    return response[ind]


def get_unique_indices(targets):
    indices = [targets.index(x) for x in set(targets)]
    return indices


def get_unique_targets(targets, *vars):
    ind_targs = get_unique_indices(targets)
    print("Number of unique flaring stars: " + str(len(ind_targs)))
    unique_vars = get_only_ind(ind_targs, *vars)

    return unique_vars


def get_only_ind(ind_targs, *vars):
    ind_vars = []
    for v in vars:
        print(len(v))
        ind_vars.append(np.array(v)[ind_targs])
    return ind_vars


def get_x_hours(total_t, X, *vars):
    ind_hours = more_than_hours(total_t, numhours=X)
    print(ind_hours)
    hours_vars = get_only_ind(ind_hours, *vars)
    return hours_vars


def more_than_hours(obs_t, numhours=20):
    obs_t = np.array(obs_t)
    return np.where(np.array(obs_t) > (numhours / 24.))[0]