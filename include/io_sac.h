#include "sacio.hpp"


// ReadSac class
class ReadSac
{
public:
    //
    ReadSac(const std::string &fname);

private:
    // necessary parameters

};

//read_sac(const std::string &fname){
//
//    Sac sac;
//    int ierr = EXIT_SUCCESS;
//
//    ierr = sac.readFile(fname);
//    if (ierr != 0) {
//        fprintf(stderr, "%s: Failed to read file\n", __func__);
//    }
//
//    // read header
//    std::string _network, _station, _channel;
//    ierr = sac.getHeader(SacHeader::String::KNETWK, _network);
//    ierr = sac.getHeader(SacHeader::String::KSTNM,  _station);
//    ierr = sac.getHeader(SacHeader::String::KCMPNM, _channel);
//    network = QString::fromStdString(_network).simplified();
//    station = QString::fromStdString(_station).simplified();
//    channel = QString::fromStdString(_channel).simplified();
//
//    ierr = sac.getHeader(SacHeader::Double::DELTA, &dt);
//    ierr = sac.getNumberOfPoints(&npts);
//
//    // read data
//    const double *data = nullptr;
//    ierr = sac.getConstantDataPtr(&data);
//    dd.resize(npts);
//    tt.resize(npts);
//    for (int i = 0; i < npts; i++){
//        dd[i] = *(data + i);
//        tt[i] = i*dt;
//    }
//
//}
