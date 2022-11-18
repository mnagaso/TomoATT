# class for storing one event data


class AttSrcRec:
    _id_src      = None
    _id_rec      = None
    _year        = None
    _month       = None
    _day         = None
    _hour        = None
    _min         = None
    _sec         = None
    _lat         = None
    _lon         = None
    _dep         = None
    _mag         = None
    _nrec        = None
    _id_event    = None
    _data_source = None
    _phase       = None
    _arr_time    = None
    _name_rec    = None

    def __init__(self,
                id_src     = None,
                id_rec     = None,
                year       = None,
                month      = None,
                day        = None,
                hour       = None,
                _min       = None,
                sec        = None,
                lat        = None,
                lon        = None,
                dep        = None,
                mag        = None,
                nrec       = None,
                id_event   = None,
                data_source= None,
                phase      = None,
                arr_time   = None,
                name_rec   = None):


        self._id_src      = id_src
        self._id_rec      = id_rec
        self._year        = year
        self._month       = month
        self._day         = day
        self._hour        = hour
        self._min         = _min
        self._sec         = sec
        self._lat         = lat
        self._lon         = lon
        self._dep         = dep
        self._mag         = mag
        self._nrec        = nrec
        self._id_event    = id_event
        self._data_source = data_source
        self._phase       = phase
        self._arr_time    = arr_time
        self._name_rec    = name_rec


def convert_to_pandas_df(event_list):
    # conveert event_list to pandas dataframe
    import pandas as pd
    import datetime

    df_ev = pd.DataFrame()

    list_id_src     = []
    list_id_rec     = []
    list_year       = []
    list_month      = []
    list_day        = []
    list_hour       = []
    list__min       = []
    list_sec        = []
    list_lat        = []
    list_lon        = []
    list_dep        = []
    list_mag        = []
    list_nrec       = []
    list_id_event   = []
    list_data_source= []
    list_phase      = []
    list_arr_time   = []
    list_name_rec   = []
    list_datetime   = []

    for ev in event_list:
        list_id_src.append(ev._id_src)
        list_id_rec.append(ev._id_rec)
        list_year.append(ev._year)
        list_month.append(ev._month)
        list_day.append(ev._day)
        list_hour.append(ev._hour)
        list__min.append(ev._min)
        list_sec.append(ev._sec)
        list_lat.append(ev._lat)
        list_lon.append(ev._lon)
        list_dep.append(ev._dep)
        list_mag.append(ev._mag)
        list_nrec.append(ev._nrec)
        list_id_event.append(ev._id_event)
        list_data_source.append(ev._data_source)
        list_phase.append(ev._phase)
        list_arr_time.append(ev._arr_time)
        list_name_rec.append(ev._name_rec)
        try:
            date_this = datetime.datetime(ev._year, ev._month, ev._day, ev._hour, ev._min, int(ev._sec))
        except:
            date_this = None

        list_datetime.append(date_this)

    # convert all the lists to pandas series
    df_ev['id_src']     = pd.Series(list_id_src)
    df_ev['id_rec']     = pd.Series(list_id_rec)
    df_ev['year']       = pd.Series(list_year)
    df_ev['month']      = pd.Series(list_month)
    df_ev['day']        = pd.Series(list_day)
    df_ev['hour']       = pd.Series(list_hour)
    df_ev['min']        = pd.Series(list__min)
    df_ev['sec']        = pd.Series(list_sec)
    df_ev['lat']        = pd.Series(list_lat)
    df_ev['lon']        = pd.Series(list_lon)
    df_ev['dep']        = pd.Series(list_dep)
    df_ev['mag']        = pd.Series(list_mag)
    df_ev['nrec']       = pd.Series(list_nrec)
    df_ev['id_event']   = pd.Series(list_id_event)
    df_ev['data_source']= pd.Series(list_data_source)
    df_ev['phase']      = pd.Series(list_phase)
    df_ev['arr_time']   = pd.Series(list_arr_time)
    df_ev['name_rec']   = pd.Series(list_name_rec)
    df_ev['datetime']   = pd.Series(list_datetime)

    return df_ev

# read file

def read_src_rec_file(fpath, two_station_names=False, data_source_flag=0, id_src_offset=0):
    #fpath = "./src_rec_test_out.dat"

    print ("read file: ", fpath)
    event_list = []
    rec_list   = []

    with open(fpath, "r") as f:
        lines = f.readlines()

        cc = 0
        nc = 0
        i_src = 0

        # parse
        for iline, line in enumerate(lines):
            #print(line)

            if line.startswith("#"):
                continue
            else:
                if cc == 0:
                    try:
                        # firstly source line is read
                        ll = line.split()

                        #src_id    = int(ll[0])
                        src_id = i_src + id_src_offset
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
                        if (nrec_tmp != 0):

                            #src = AttEvent(src_id, src_year, src_month, src_day, src_hour, src_min, src_sec, src_lat, src_lon, src_dep, src_mag, src_nrec, src_id_event, data_source_flag)

                            src = AttSrcRec(src_id, None, src_year, src_month, src_day, src_hour, src_min, src_sec, src_lat, src_lon, src_dep, src_mag, src_nrec, src_id_event, data_source_flag, None, None, None, None)
                            event_list.append(src)

                            cc+=1
                        else:
                            pass
                    except:
                        pass
                else:
                    try:
                        # read rec line
                        ll = line.split()

                        if(not two_station_names):
                            #src_id   = int(ll[0])
                            src_id = i_src + id_src_offset
                            rec_id   = int(ll[1])
                            rec_name = ll[2]
                            rec_lat  = float(ll[3])
                            rec_lon  = float(ll[4])
                            rec_elev  = float(ll[5])
                            rec_phase = ll[6]
                            rec_arr_time = float(ll[8])
                        else:
                            #src_id   = int(ll[0])
                            src_id = i_src + id_src_offset
                            rec_id   = int(ll[1])
                            rec_name = ll[2] + "_" + ll[3]
                            rec_lat  = float(ll[4])
                            rec_lon  = float(ll[5])
                            rec_elev  = float(ll[6])
                            rec_phase = ll[7]
                            rec_arr_time = float(ll[9])

                        # store rec
                        #rec = AttArrival(src_id, rec_id, rec_name, rec_lat, rec_lon, rec_elev, rec_phase, rec_epi_dist, rec_arr_time)
                        #event_list[i_src].add_rec(rec)
                        rec = AttSrcRec(src_id, rec_id, None, None, None, None, None, None, rec_lat, rec_lon, rec_elev, None, None, None, data_source_flag, rec_phase, rec_arr_time, rec_name)
                        rec_list.append(rec)

                        nc+=1
                    except:
                        print("error in line: ", iline)
                        print("error in line: " + line)
                        #return None

                    cc+=1

                    if cc > nrec_tmp:
                        cc = 0
                        i_src += 1

                        if nc == 0:
                            print("error: no rec found")
                            # erase last event
                            event_list.pop()
                            i_src -= 1

                        nc = 0


    # return length of event_list
    print("number of events: ", len(event_list))
    print("number of recs: ", len(rec_list))

    df_ev = convert_to_pandas_df(event_list)
    df_rec = convert_to_pandas_df(rec_list)


    return df_ev, df_rec


if __name__ == "__main__":
    event_list = read_src_rec_file("./src_rec_test_out.dat")
    print(event_list[0].rec_list[0].name_rec)
    print(event_list[0].rec_list[0].arr_time)
    print(event_list[0].rec_list[0].id_rec)
    print(event_list[0].rec_list[0].id_src)
    print(event_list[0].rec_list[0].lat)
    print(event_list[0].rec_list[0].lon)
    print(event_list[0].rec_list[0].dep)
    print(event_list[0].rec_list[0].phase)
    print(event_list[0].rec_list[0].arr_time)
    print(event_list[0].rec_list[1].name_rec)
    print(event_list[0].rec_list[1].arr_time)
    print(event_list[0].rec_list[1].id_rec)
    print(event_list[0].rec_list[1].id_src)
    print(event_list[0].rec_list[1].lat)
    print(event_list[0].rec_list[1].lon)
    print(event_list[0].rec_list[1].dep)
    print(event_list[0].rec_list[1].phase)
    print(event_list[0].rec_list[1].arr_time)
    print(event_list[0].rec_list[2].name_rec)