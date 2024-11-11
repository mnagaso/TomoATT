# %%
# Initialization, class definitions, and declarations
from obspy import UTCDateTime

class Event():      # Class of earthquake
    def __init__(self):
        self.name = "nan"       # Earthquake name, recommended as "earthquake"
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.dep = 0.0
        self.mag = 0.0
        self.ortime = UTCDateTime(1999,1,1,0,0,0)
        self.Nt = 0             # Number of the absolute traveltime of earthquake
        self.Ncs_dt = 0         # Number of the common source differential traveltime of earthquake
        self.Ncr_dt = 0         # Number of the common receiver differential traveltime of earthquake
        self.t = {}             # stname1+phase -> (stname1, phase, time, data_weight)
        self.cs_dt = {}         # stname1 + stname2 + phase -> (stname1, stname2, phase, dif_time, data_weight)
        self.cr_dt = {}         # stname1 + evname2 + phase -> (stname1, evname2, phase, dif_time, data_weight)
        self.azi_gap = 360.0    # The max azimuthal gap of this earthquake
        self.misfit = {}        # Traveltime residual of the data, the difference between real data and synthetic data, used for evaluation. stname or stname1+stname2 or stname1+evname2 -> residual
        self.tag    = {}        # Additional tags for the earthquake, e.g., azi_gap, weight. (azimuthal gap, weight of the earthquake)

class Station():
    def __init__(self):
        self.name = "nan"       # Station name, recommended as network.station_name
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.ele = 0.0
        self.tag = {}           # Additional tags for the station, e.g., weight, used to describe the attributes of the station, such as station weight: weight

# %% [markdown]
# Functions: some basic auxiliary functions for processing data

# %%
# Function: cal_dis(lat1, lon1, lat2, lon2) (km), cal_azimuth(lat1, lon1, lat2, lon2) (°) compute epicentral distance (km) and azimuth (degree)

def cal_dis(lat1, lon1, lat2, lon2, R = 6371):
    latitude1 = (math.pi/180)*lat1
    latitude2 = (math.pi/180)*lat2
    longitude1 = (math.pi/180)*lon1
    longitude2 = (math.pi/180)*lon2
    # Therefore, the spherical distance between points A and B is: {arccos[sinb*siny+cosb*cosy*cos(a-x)]}*R
    # Earth radius
    if((lat1-lat2)**2+(lon1-lon2)**2<0.000001):
        return 0

    d = math.acos(math.sin(latitude1)*math.sin(latitude2)+ math.cos(latitude1)*math.cos(latitude2)*math.cos(longitude2-longitude1))/math.pi*180
    return d * 2 * math.pi * R / 360

def cal_azimuth(lat1, lon1, lat2, lon2):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180

    y = math.sin(lon2_rad - lon1_rad) * math.cos(lat2_rad)
    x = math.cos(lat1_rad) * math.sin(lat2_rad) - math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(lon2_rad - lon1_rad)
    brng = math.atan2(y, x) * 180 / math.pi
    if((lat1-lat2)**2+(lon1-lon2)**2<0.0001):
        return 0
    return float((brng + 360.0) % 360.0)

# %%
# Function: Coordinate rotation rotate_src_rec(ev_info, st_info, theta0, phi0, psi): rotate to the new coordinate, satisfying the center r0, t0, p0 -> r0, 0, 0 and an anticlockwise angle psi

RAD2DEG = 180/np.pi
DEG2RAD = np.pi/180
R_earth = 6371.0

# Spherical coordinates to Cartesian coordinate
def rtp2xyz(r, theta, phi):
    x = r * np.cos(theta*DEG2RAD) * np.cos(phi*DEG2RAD)
    y = r * np.cos(theta*DEG2RAD) * np.sin(phi*DEG2RAD)
    z = r * np.sin(theta*DEG2RAD)
    return (x, y, z)

# Cartesian coordinates to Spherical coordinate
def xyz2rtp(x, y, z):
    # theta: -90~90;  phi: -180~180
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arcsin(z/r)
    phi = np.arcsin(y/r/np.cos(theta))

    idx = np.where((phi > 0) & (x*y < 0))
    phi[idx] = np.pi - phi[idx]
    idx = np.where((phi < 0) & (x*y > 0))
    phi[idx] = -np.pi - phi[idx]

    return (r, theta*RAD2DEG, phi*RAD2DEG)

# Anti-clockwise rotation along x-axis
def rotate_x(x, y, z, theta):
    new_x = x
    new_y = y *  np.cos(theta*DEG2RAD) + z * -np.sin(theta*DEG2RAD)
    new_z = y *  np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
    return (new_x, new_y, new_z)

# Anti-clockwise rotation along y-axis
def rotate_y(x, y, z, theta):
    new_x = x *  np.cos(theta*DEG2RAD) + z *  np.sin(theta*DEG2RAD)
    new_y = y
    new_z = x * -np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
    return (new_x, new_y, new_z)

# Anti-clockwise rotation along z-axis
def rotate_z(x, y, z, theta):
    new_x = x *  np.cos(theta*DEG2RAD) + y * -np.sin(theta*DEG2RAD)
    new_y = x *  np.sin(theta*DEG2RAD) + y *  np.cos(theta*DEG2RAD)
    new_z = z
    return (new_x, new_y, new_z)

# Spherical Rotation

# Rotate to the new coordinate, satisfying the center r0, t0, p0 -> r0, 0, 0 and an anticlockwise angle psi
def rtp_rotation(t, p, theta0, phi0, psi):
    # Step 1: r, t, p -> x, y, z
    (x, y, z) = rtp2xyz(1.0, t, p)

    # Step 2: Anti-clockwise rotation with -phi0 along z-axis: r0, t0, p0 -> r0, t0, 0
    (x, y, z) = rotate_z(x, y, z, -phi0)

    # Step 3: Anti-clockwise rotation with theta0 along y-axis: r0, t0, 0 -> r0, 0, 0
    (x, y, z) = rotate_y(x, y, z, theta0)

    # Step 4: Anti-clockwise rotation with psi along x-axis
    (x, y, z) = rotate_x(x, y, z, psi)

    # Step 5: x, y, z -> r, t, p
    (new_r, new_t, new_p) = xyz2rtp(x, y, z)

    return (new_t, new_p)

def rtp_rotation_reverse(new_t, new_p, theta0, phi0, psi):
    # Step 1: r, t, p -> x, y, z
    (x, y, z) = rtp2xyz(1.0, new_t, new_p)

    # Step 2: Anti-clockwise rotation with -psi along x-axis
    (x, y, z) = rotate_x(x, y, z, -psi)

    # Step 3: Anti-clockwise rotation with -theta0 along y-axis: r0, 0, 0 -> r0, t0, 0
    (x, y, z) = rotate_y(x, y, z, -theta0)

    # Step 4: Anti-clockwise rotation with phi0 along z-axis: r0, t0, 0 -> r0, t0, p0
    (x, y, z) = rotate_z(x, y, z, phi0)

    # Step 5: x, y, z -> r, t, p
    (r, t, p) = xyz2rtp(x, y, z)

    return (t, p)

def rotate_src_rec(ev_info, st_info, theta0, phi0, psi):
    ev_info_rotate = {}
    st_info_rotate = {}

    # Rotate earthquakes
    for key_ev in ev_info:
        ev = ev_info[key_ev]
        ev_lat = np.array([ev.lat])
        ev_lon = np.array([ev.lon])
        (ev_lat, ev_lon,) = rtp_rotation(ev_lat, ev_lon, theta0, phi0, psi)
        ev.lat = ev_lat[0]
        ev.lon = ev_lon[0]
        ev_info_rotate[key_ev] = ev

    # Rotate stations
    for key_st in st_info:
        st = st_info[key_st]
        st_lat = np.array([st.lat])
        st_lon = np.array([st.lon])
        (st_lat, st_lon) = rtp_rotation(st_lat, st_lon, theta0, phi0, psi)
        st.lat = st_lat[0]
        st.lon = st_lon[0]
        st_info_rotate[key_st] = st

    return (ev_info_rotate, st_info_rotate)

def rotate_src_rec_reverse(ev_info_rotate, st_info_rotate, theta0, phi0, psi):
    ev_info = {}
    st_info = {}

    # Rotate earthquakes
    for key_ev in ev_info_rotate:
        ev = ev_info_rotate[key_ev]
        ev_lat = np.array([ev.lat])
        ev_lon = np.array([ev.lon])
        (ev_lat, ev_lon,) = rtp_rotation_reverse(ev_lat, ev_lon, theta0, phi0, psi)
        ev.lat = ev_lat[0]
        ev.lon = ev_lon[0]
        ev_info[key_ev] = ev

    # Rotate stations
    for key_st in st_info_rotate:
        st = st_info_rotate[key_st]
        st_lat = np.array([st.lat])
        st_lon = np.array([st.lon])
        (st_lat, st_lon) = rtp_rotation_reverse(st_lat, st_lon, theta0, phi0, psi)
        st.lat = st_lat[0]
        st.lon = st_lon[0]
        st_info[key_st] = st

    return (ev_info, st_info)

# %%
# Function: linear_regression(X, Y) compute linear regression
def linear_regression(X, Y):
    slope, intercept = np.polyfit(X, Y, deg=1)
    fitted_values = slope * X + intercept
    residual = Y - fitted_values
    SEE = np.std(residual)
    return (slope, intercept, SEE)

# %% [markdown]
# Functions: obtain target information from ev_info and st_info

# %%
# Function: data_lon_lat_dep_wt_ev(ev_info) output the [lon, lat, dep, weight] of the earthquake
def data_lon_lat_dep_wt_ev(ev_info):
    lat = []
    lon = []
    dep = []
    weight = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        try:
            weight.append(ev_info[key].tag["weight"])
        except:
            weight.append(1.0)
    return [np.array(lon), np.array(lat), np.array(dep), np.array(weight)]

# %%
# Function: data_ev_loc(ev_info) output the [lon, lat, dep, ortime] of the earthquake
def data_ev_loc(ev_info):
    lat = []
    lon = []
    dep = []
    ortime = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        ortime.append(ev_info[key].ortime.timestamp)
    return [np.array(lon), np.array(lat), np.array(dep), np.array(ortime)]

# %%
# Function: data_lon_lat_ele_wt_st(ev_info, st_info) output the [lon, lat, dep, weight] of the station
def data_lon_lat_ele_wt_st(ev_info, st_info):
    names = {}
    lat = []
    lon = []
    ele = []
    weight  = []
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            names[ev_info[key_ev].t[key_t][0]] = 1

        for key_t in ev_info[key_ev].cs_dt:
            names[ev_info[key_ev].cs_dt[key_t][0]] = 1
            names[ev_info[key_ev].cs_dt[key_t][1]] = 1

        for key_t in ev_info[key_ev].cr_dt:
            names[ev_info[key_ev].cr_dt[key_t][0]] = 1

    for name in names:  # Only output the station which has data
        lat.append(st_info[name].lat)
        lon.append(st_info[name].lon)
        ele.append(st_info[name].ele)
        try:
            weight.append(st_info[name].tag["weight"])
        except:
            weight.append(1.0)
    return [np.array(lon), np.array(lat), np.array(ele), np.array(weight)]

# %%
# Function: data_dis_time(ev_info) output the [dis, time] of all data
def data_dis_time(ev_info, st_info):
    all_dis = []
    all_time = []
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
            dep_st = st_info[ev_info[key_ev].t[key_t][0]].ele
            dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            all_dis.append(dis)
            all_time.append(ev_info[key_ev].t[key_t][2])

    return [np.array(all_dis), np.array(all_time)]

# %%
# Function: data_epidis_time(ev_info) output the [epidis, time] of all data
def data_epidis_time(ev_info, st_info):
    all_dis = []
    all_time = []
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
            dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            all_dis.append(dis)
            all_time.append(ev_info[key_ev].t[key_t][2])

    return [np.array(all_dis), np.array(all_time)]

# %%
# Function: data_cs_dt(ev_info) output the [cs_dt] of all data
def data_cs_dt(ev_info):
    all_time = []
    for key_ev in ev_info:
        for key_dt in ev_info[key_ev].cs_dt:
            all_time.append(ev_info[key_ev].cs_dt[key_dt][3])

    return np.array(all_time)

# %%
# Function: data_dis_time_phase(ev_info, st_info, phase_list) output (dis, time) of each given phase
def data_dis_time_phase(ev_info, st_info, phase_list):
    all_dis  = {}
    all_time = {}
    for phase in phase_list:
        all_dis[phase] = []
        all_time[phase] = []

    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        for key_t in ev_info[key_ev].t:
            if ev_info[key_ev].t[key_t][1] in phase_list:
                lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
                lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
                dep_st = st_info[ev_info[key_ev].t[key_t][0]].ele
                dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
                all_dis[ev_info[key_ev].t[key_t][1]].append(dis)
                all_time[ev_info[key_ev].t[key_t][1]].append(ev_info[key_ev].t[key_t][2])

    for phase in phase_list:
        all_dis[phase] = np.array(all_dis[phase])
        all_time[phase] = np.array(all_time[phase])

    return [all_dis, all_time]

# %%
# Function: data_line(ev_info, st_info) output the lines connecting station and earthquake
def data_line(ev_info, st_info):
    line_x = []
    line_y = []

    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t# %%
# Initialization, class definitions, and declarations
from obspy import UTCDateTime

class Event():      # Class of earthquake
    def __init__(self):
        self.name = "nan"       # Earthquake name, recommended as "earthquake"
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.dep = 0.0
        self.mag = 0.0
        self.ortime = UTCDateTime(1999,1,1,0,0,0)
        self.Nt = 0             # Number of the absolute traveltime of earthquake
        self.Ncs_dt = 0         # Number of the common source differential traveltime of earthquake
        self.Ncr_dt = 0         # Number of the common receiver differential traveltime of earthquake
        self.t = {}             # stname1+phase -> (stname1, phase, time, data_weight)
        self.cs_dt = {}         # stname1 + stname2 + phase -> (stname1, stname2, phase, dif_time, data_weight)
        self.cr_dt = {}         # stname1 + evname2 + phase -> (stname1, evname2, phase, dif_time, data_weight)
        self.azi_gap = 360.0    # The max azimuthal gap of this earthquake
        self.misfit = {}        # Traveltime residual of the data, the difference between real data and synthetic data, used for evaluation. stname or stname1+stname2 or stname1+evname2 -> residual
        self.tag    = {}        # Additional tags for the earthquake, e.g., azi_gap, weight. (azimuthal gap, weight of the earthquake)

class Station():
    def __init__(self):
        self.name = "nan"       # Station name, recommended as network.station_name
        self.id = -1
        self.lat = 0.0
        self.lon = 0.0
        self.ele = 0.0
        self.tag = {}           # Additional tags for the station, e.g., weight, used to describe the attributes of the station, such as station weight: weight

# %% [markdown]
# Functions: some basic auxiliary functions for processing data

# %%
# Function: cal_dis(lat1, lon1, lat2, lon2) (km), cal_azimuth(lat1, lon1, lat2, lon2) (°) compute epicentral distance (km) and azimuth (degree)

def cal_dis(lat1, lon1, lat2, lon2, R = 6371):
    latitude1 = (math.pi/180)*lat1
    latitude2 = (math.pi/180)*lat2
    longitude1 = (math.pi/180)*lon1
    longitude2 = (math.pi/180)*lon2
    # Therefore, the spherical distance between points A and B is: {arccos[sinb*siny+cosb*cosy*cos(a-x)]}*R
    # Earth radius
    if((lat1-lat2)**2+(lon1-lon2)**2<0.000001):
        return 0

    d = math.acos(math.sin(latitude1)*math.sin(latitude2)+ math.cos(latitude1)*math.cos(latitude2)*math.cos(longitude2-longitude1))/math.pi*180
    return d * 2 * math.pi * R / 360

def cal_azimuth(lat1, lon1, lat2, lon2):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180

    y = math.sin(lon2_rad - lon1_rad) * math.cos(lat2_rad)
    x = math.cos(lat1_rad) * math.sin(lat2_rad) - math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(lon2_rad - lon1_rad)
    brng = math.atan2(y, x) * 180 / math.pi
    if((lat1-lat2)**2+(lon1-lon2)**2<0.0001):
        return 0
    return float((brng + 360.0) % 360.0)

# %%
# Function: Coordinate rotation rotate_src_rec(ev_info, st_info, theta0, phi0, psi): rotate to the new coordinate, satisfying the center r0, t0, p0 -> r0, 0, 0 and an anticlockwise angle psi

RAD2DEG = 180/np.pi
DEG2RAD = np.pi/180
R_earth = 6371.0

# Spherical coordinates to Cartesian coordinate
def rtp2xyz(r, theta, phi):
    x = r * np.cos(theta*DEG2RAD) * np.cos(phi*DEG2RAD)
    y = r * np.cos(theta*DEG2RAD) * np.sin(phi*DEG2RAD)
    z = r * np.sin(theta*DEG2RAD)
    return (x, y, z)

# Cartesian coordinates to Spherical coordinate
def xyz2rtp(x, y, z):
    # theta: -90~90;  phi: -180~180
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arcsin(z/r)
    phi = np.arcsin(y/r/np.cos(theta))

    idx = np.where((phi > 0) & (x*y < 0))
    phi[idx] = np.pi - phi[idx]
    idx = np.where((phi < 0) & (x*y > 0))
    phi[idx] = -np.pi - phi[idx]

    return (r, theta*RAD2DEG, phi*RAD2DEG)

# Anti-clockwise rotation along x-axis
def rotate_x(x, y, z, theta):
    new_x = x
    new_y = y *  np.cos(theta*DEG2RAD) + z * -np.sin(theta*DEG2RAD)
    new_z = y *  np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
    return (new_x, new_y, new_z)

# Anti-clockwise rotation along y-axis
def rotate_y(x, y, z, theta):
    new_x = x *  np.cos(theta*DEG2RAD) + z *  np.sin(theta*DEG2RAD)
    new_y = y
    new_z = x * -np.sin(theta*DEG2RAD) + z *  np.cos(theta*DEG2RAD)
    return (new_x, new_y, new_z)

# Anti-clockwise rotation along z-axis
def rotate_z(x, y, z, theta):
    new_x = x *  np.cos(theta*DEG2RAD) + y * -np.sin(theta*DEG2RAD)
    new_y = x *  np.sin(theta*DEG2RAD) + y *  np.cos(theta*DEG2RAD)
    new_z = z
    return (new_x, new_y, new_z)

# Spherical Rotation

# Rotate to the new coordinate, satisfying the center r0, t0, p0 -> r0, 0, 0 and an anticlockwise angle psi
def rtp_rotation(t, p, theta0, phi0, psi):
    # Step 1: r, t, p -> x, y, z
    (x, y, z) = rtp2xyz(1.0, t, p)

    # Step 2: Anti-clockwise rotation with -phi0 along z-axis: r0, t0, p0 -> r0, t0, 0
    (x, y, z) = rotate_z(x, y, z, -phi0)

    # Step 3: Anti-clockwise rotation with theta0 along y-axis: r0, t0, 0 -> r0, 0, 0
    (x, y, z) = rotate_y(x, y, z, theta0)

    # Step 4: Anti-clockwise rotation with psi along x-axis
    (x, y, z) = rotate_x(x, y, z, psi)

    # Step 5: x, y, z -> r, t, p
    (new_r, new_t, new_p) = xyz2rtp(x, y, z)

    return (new_t, new_p)

def rtp_rotation_reverse(new_t, new_p, theta0, phi0, psi):
    # Step 1: r, t, p -> x, y, z
    (x, y, z) = rtp2xyz(1.0, new_t, new_p)

    # Step 2: Anti-clockwise rotation with -psi along x-axis
    (x, y, z) = rotate_x(x, y, z, -psi)

    # Step 3: Anti-clockwise rotation with -theta0 along y-axis: r0, 0, 0 -> r0, t0, 0
    (x, y, z) = rotate_y(x, y, z, -theta0)

    # Step 4: Anti-clockwise rotation with phi0 along z-axis: r0, t0, 0 -> r0, t0, p0
    (x, y, z) = rotate_z(x, y, z, phi0)

    # Step 5: x, y, z -> r, t, p
    (r, t, p) = xyz2rtp(x, y, z)

    return (t, p)

def rotate_src_rec(ev_info, st_info, theta0, phi0, psi):
    ev_info_rotate = {}
    st_info_rotate = {}

    # Rotate earthquakes
    for key_ev in ev_info:
        ev = ev_info[key_ev]
        ev_lat = np.array([ev.lat])
        ev_lon = np.array([ev.lon])
        (ev_lat, ev_lon,) = rtp_rotation(ev_lat, ev_lon, theta0, phi0, psi)
        ev.lat = ev_lat[0]
        ev.lon = ev_lon[0]
        ev_info_rotate[key_ev] = ev

    # Rotate stations
    for key_st in st_info:
        st = st_info[key_st]
        st_lat = np.array([st.lat])
        st_lon = np.array([st.lon])
        (st_lat, st_lon) = rtp_rotation(st_lat, st_lon, theta0, phi0, psi)
        st.lat = st_lat[0]
        st.lon = st_lon[0]
        st_info_rotate[key_st] = st

    return (ev_info_rotate, st_info_rotate)

def rotate_src_rec_reverse(ev_info_rotate, st_info_rotate, theta0, phi0, psi):
    ev_info = {}
    st_info = {}

    # Rotate earthquakes
    for key_ev in ev_info_rotate:
        ev = ev_info_rotate[key_ev]
        ev_lat = np.array([ev.lat])
        ev_lon = np.array([ev.lon])
        (ev_lat, ev_lon,) = rtp_rotation_reverse(ev_lat, ev_lon, theta0, phi0, psi)
        ev.lat = ev_lat[0]
        ev.lon = ev_lon[0]
        ev_info[key_ev] = ev

    # Rotate stations
    for key_st in st_info_rotate:
        st = st_info_rotate[key_st]
        st_lat = np.array([st.lat])
        st_lon = np.array([st.lon])
        (st_lat, st_lon) = rtp_rotation_reverse(st_lat, st_lon, theta0, phi0, psi)
        st.lat = st_lat[0]
        st.lon = st_lon[0]
        st_info[key_st] = st

    return (ev_info, st_info)

# %%
# Function: linear_regression(X, Y) compute linear regression
def linear_regression(X, Y):
    slope, intercept = np.polyfit(X, Y, deg=1)
    fitted_values = slope * X + intercept
    residual = Y - fitted_values
    SEE = np.std(residual)
    return (slope, intercept, SEE)

# %% [markdown]
# Functions: obtain target information from ev_info and st_info

# %%
# Function: data_lon_lat_dep_wt_ev(ev_info) output the [lon, lat, dep, weight] of the earthquake
def data_lon_lat_dep_wt_ev(ev_info):
    lat = []
    lon = []
    dep = []
    weight = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        try:
            weight.append(ev_info[key].tag["weight"])
        except:
            weight.append(1.0)
    return [np.array(lon), np.array(lat), np.array(dep), np.array(weight)]

# %%
# Function: data_ev_loc(ev_info) output the [lon, lat, dep, ortime] of the earthquake
def data_ev_loc(ev_info):
    lat = []
    lon = []
    dep = []
    ortime = []
    for key in ev_info:
        lat.append(ev_info[key].lat)
        lon.append(ev_info[key].lon)
        dep.append(ev_info[key].dep)
        ortime.append(ev_info[key].ortime.timestamp)
    return [np.array(lon), np.array(lat), np.array(dep), np.array(ortime)]

# %%
# Function: data_lon_lat_ele_wt_st(ev_info, st_info) output the [lon, lat, dep, weight] of the station
def data_lon_lat_ele_wt_st(ev_info, st_info):
    names = {}
    lat = []
    lon = []
    ele = []
    weight  = []
    for key_ev in ev_info:
        for key_t in ev_info[key_ev].t:
            names[ev_info[key_ev].t[key_t][0]] = 1

        for key_t in ev_info[key_ev].cs_dt:
            names[ev_info[key_ev].cs_dt[key_t][0]] = 1
            names[ev_info[key_ev].cs_dt[key_t][1]] = 1

        for key_t in ev_info[key_ev].cr_dt:
            names[ev_info[key_ev].cr_dt[key_t][0]] = 1

    for name in names:  # Only output the station which has data
        lat.append(st_info[name].lat)
        lon.append(st_info[name].lon)
        ele.append(st_info[name].ele)
        try:
            weight.append(st_info[name].tag["weight"])
        except:
            weight.append(1.0)
    return [np.array(lon), np.array(lat), np.array(ele), np.array(weight)]

# %%
# Function: data_dis_time(ev_info) output the [dis, time] of all data
def data_dis_time(ev_info, st_info):
    all_dis = []
    all_time = []
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
            dep_st = st_info[ev_info[key_ev].t[key_t][0]].ele
            dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            all_dis.append(dis)
            all_time.append(ev_info[key_ev].t[key_t][2])

    return [np.array(all_dis), np.array(all_time)]

# %%
# Function: data_epidis_time(ev_info) output the [epidis, time] of all data
def data_epidis_time(ev_info, st_info):
    all_dis = []
    all_time = []
    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
            lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
            dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
            all_dis.append(dis)
            all_time.append(ev_info[key_ev].t[key_t][2])

    return [np.array(all_dis), np.array(all_time)]

# %%
# Function: data_cs_dt(ev_info) output the [cs_dt] of all data
def data_cs_dt(ev_info):
    all_time = []
    for key_ev in ev_info:
        for key_dt in ev_info[key_ev].cs_dt:
            all_time.append(ev_info[key_ev].cs_dt[key_dt][3])

    return np.array(all_time)

# %%
# Function: data_dis_time_phase(ev_info, st_info, phase_list) output (dis, time) of each given phase
def data_dis_time_phase(ev_info, st_info, phase_list):
    all_dis  = {}
    all_time = {}
    for phase in phase_list:
        all_dis[phase] = []
        all_time[phase] = []

    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        dep_ev = ev_info[key_ev].dep
        for key_t in ev_info[key_ev].t:
            if ev_info[key_ev].t[key_t][1] in phase_list:
                lat_st = st_info[ev_info[key_ev].t[key_t][0]].lat
                lon_st = st_info[ev_info[key_ev].t[key_t][0]].lon
                dep_st = st_info[ev_info[key_ev].t[key_t][0]].ele
                dis = cal_dis(lat_ev, lon_ev, lat_st, lon_st)
                all_dis[ev_info[key_ev].t[key_t][1]].append(dis)
                all_time[ev_info[key_ev].t[key_t][1]].append(ev_info[key_ev].t[key_t][2])

    for phase in phase_list:
        all_dis[phase] = np.array(all_dis[phase])
        all_time[phase] = np.array(all_time[phase])

    return [all_dis, all_time]

# %%
# Function: data_line(ev_info, st_info) output the lines connecting station and earthquake
def data_line(ev_info, st_info):
    line_x = []
    line_y = []

    for key_ev in ev_info:
        lat_ev = ev_info[key_ev].lat
        lon_ev = ev_info[key_ev].lon
        for key_t in ev_info[key_ev].t:
            lat_st = st_info[ev_info[key_ev].t[key_t