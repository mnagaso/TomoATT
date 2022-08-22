# class for storing one event data

class AttArrival:
    def __init__(self, id_src,id_rec,name_rec,lat,lon,dep,phase,epi_dist,arr_time):
        self.id_src = id_src
        self.id_rec = id_rec
        self.name_rec = name_rec
        self.lat = lat
        self.lon = lon
        self.dep = dep
        self.phase = phase
        self.epi_dist = epi_dist
        self.arr_time = arr_time

class AttEvent:
    def __init__(self, _id, year, month, day, hour, _min, sec, lat, lon, dep, mag, nrec, id_event):
        self.id = _id
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.min = _min
        self.sec = sec
        self.lat = lat
        self.lon = lon
        self.dep = dep
        self.mag = mag
        self.nrec = nrec
        self.id_event = id_event
        self.rec_list = []

    def add_rec(self, rec):
        self.rec_list.append(rec)

# read file

def read_src_rec_file(fpath):
    #fpath = "./src_rec_test_out.dat"

    event_list = []

    with open(fpath, "r") as f:
        lines = f.readlines()

        cc = 0
        i_src = 0

        # parse
        for line in lines:
            if line.startswith("#"):
                continue
            else:
                if cc == 0:
                    # firstly source line is read
                    ll = line.split()
                    src_id    = int(ll[0])
                    src_year  = int(ll[1])
                    src_month = int(ll[2])
                    src_day   = int(ll[3])
                    src_hour  = int(ll[4])
                    src_min   = int(ll[5])
                    src_sec   = float(ll[6])
                    src_lat   = float(ll[7])
                    src_lon   = float(ll[8])
                    src_dep   = float(ll[9])
                    src_mag   = float(ll[10])
                    src_nrec  = int(ll[11])
                    src_id_event = ll[12]

                    nrec_tmp = src_nrec

                    # store source
                    if nrec_tmp != 0:
                        src = AttEvent(src_id, src_year, src_month, src_day, src_hour, src_min, src_sec, src_lat, src_lon, src_dep, src_mag, src_nrec, src_id_event)
                        event_list.append(src)
                        cc+=1
                    else:
                        pass
                else:
                    # read rec line
                    ll = line.split()
                    src_id   = int(ll[0])
                    rec_id   = int(ll[1])
                    rec_name = ll[2]
                    rec_lat  = float(ll[3])
                    rec_lon  = float(ll[4])
                    rec_dep  = float(ll[5])
                    rec_phase = ll[6]
                    rec_epi_dist = float(ll[7])
                    rec_arr_time = float(ll[8])

                    # store rec
                    rec = AttArrival(src_id, rec_id, rec_name, rec_lat, rec_lon, rec_dep, rec_phase, rec_epi_dist, rec_arr_time)
                    event_list[i_src].add_rec(rec)

                    cc+=1

                    if cc > nrec_tmp:
                        cc = 0
                        i_src += 1


    # return length of event_list
    print("number of events: ", len(event_list))

    return event_list


if __name__ is "__main__":
    event_list = read_src_rec_file("./src_rec_test_out.dat")
    print(event_list[0].rec_list[0].name_rec)
    print(event_list[0].rec_list[0].epi_dist)
    print(event_list[0].rec_list[0].arr_time)
    print(event_list[0].rec_list[0].id_rec)
    print(event_list[0].rec_list[0].id_src)
    print(event_list[0].rec_list[0].lat)
    print(event_list[0].rec_list[0].lon)
    print(event_list[0].rec_list[0].dep)
    print(event_list[0].rec_list[0].phase)
    print(event_list[0].rec_list[0].epi_dist)
    print(event_list[0].rec_list[0].arr_time)
    print(event_list[0].rec_list[1].name_rec)
    print(event_list[0].rec_list[1].epi_dist)
    print(event_list[0].rec_list[1].arr_time)
    print(event_list[0].rec_list[1].id_rec)
    print(event_list[0].rec_list[1].id_src)
    print(event_list[0].rec_list[1].lat)
    print(event_list[0].rec_list[1].lon)
    print(event_list[0].rec_list[1].dep)
    print(event_list[0].rec_list[1].phase)
    print(event_list[0].rec_list[1].epi_dist)
    print(event_list[0].rec_list[1].arr_time)
    print(event_list[0].rec_list[2].name_rec)
    print(event_list[0].rec_list[2].epi_dist)