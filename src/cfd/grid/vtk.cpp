#include "vtk.hpp"
#include <filesystem>
#include <iomanip>
#include <sstream>

using namespace cfd;

namespace fs = std::filesystem;

namespace {

template<typename C> void set_ndata(const std::vector<C>& data, size_t& ndata) {
    if (ndata == 0) {
        ndata = data.size();
    }
    if (ndata > data.size()) {
        throw std::runtime_error("data length is invalid");
    }
}

std::vector<size_t> find_strings_in_stream(std::string fname, const std::vector<std::string>& str) {
    std::ifstream fs(fname);
    std::string out(std::istreambuf_iterator<char>{fs}, std::istreambuf_iterator<char>{});
    std::vector<size_t> ret;
    for (size_t i = 0; i < str.size(); ++i)
        ret.push_back(out.find(str[i]));
    return ret;
}

// std::ostream& vtk_string(std::ostream& os, const Point& p) {
//     os << p[0] << " " << p[1] << " " << p[2];
//     return os;
// }
//
// std::ostream& vtk_string(std::ostream& s, double v) {
//     s << v;
//     return s;
// }
//
// std::ostream& vtk_string(std::ostream& s, const std::vector<double>& v) {
//     for (double d : v) {
//         s << d << std::endl;
//     }
//     return s;
// }
//
// std::ostream& vtk_string(std::ostream& s, const std::vector<Vector>& v) {
//     for (Vector d : v) {
//         s << d[0] << " " << d[1] << " " << d[2] << std::endl;
//     }
//     return s;
// }

std::ostream& vtk_string(std::ostream& s, std::vector<double>::const_iterator v, size_t ndata) {
    for (size_t i = 0; i < ndata; ++i) {
        s << *v++ << std::endl;
    }
    return s;
}

std::ostream& vtk_string(std::ostream& s, std::vector<double>::const_iterator x, std::vector<double>::const_iterator y,
                         size_t ndata) {
    for (size_t i = 0; i < ndata; ++i) {
        s << *x++ << " " << *y++ << " 0" << std::endl;
    }
    return s;
}

std::ostream& vtk_string(std::ostream& s, std::vector<double>::const_iterator x, std::vector<double>::const_iterator y,
                         std::vector<double>::const_iterator z, size_t ndata) {
    for (size_t i = 0; i < ndata; ++i) {
        s << *x++ << " " << *y++ << " " << *z++ << std::endl;
    }
    return s;
}

std::ostream& vtk_string(std::ostream& s, std::vector<Vector>::const_iterator v, size_t ndata) {
    for (size_t i = 0; i < ndata; ++i) {
        s << (*v).x << " " << (*v).y << " " << (*v).z << std::endl;
        v++;
    }
    return s;
}

void create_directory(std::string path, bool purge) {
    if (fs::is_directory(path)) {
        if (purge) {
            fs::remove_all(path);
        } else {
            return;
        }
    }
    fs::create_directory(path);
}

void add_data_stream(std::vector<double>::const_iterator begin, size_t ndata, std::string data_cap, std::ostream& fs) {
    fs << "SCALARS " << data_cap << " double 1" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    vtk_string(fs, begin, ndata);
}

void add_vector_stream(std::vector<Vector>::const_iterator begin, size_t ndata, std::string data_cap,
                       std::ostream& fs) {
    fs << "SCALARS " << data_cap << " double 3" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    vtk_string(fs, begin, ndata);
}

void add_vector_stream(std::vector<double>::const_iterator x, std::vector<double>::const_iterator y, size_t ndata,
                       std::string data_cap, std::ostream& fs) {
    fs << "SCALARS " << data_cap << " double 3" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    vtk_string(fs, x, y, ndata);
}

void add_vector_stream(std::vector<double>::const_iterator x, std::vector<double>::const_iterator y,
                       std::vector<double>::const_iterator z, size_t ndata, std::string data_cap, std::ostream& fs) {

    fs << "SCALARS " << data_cap << " double 3" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    vtk_string(fs, x, y, z, ndata);
}

void append_cell_data_header(size_t data_size, std::ostream& os) {
    os << "CELL_DATA " << data_size << std::endl;
}

void append_point_data_header(size_t data_size, std::ostream& os) {
    os << "POINT_DATA " << data_size << std::endl;
}

void write_cell_string(std::string fname, const std::string& data, size_t ndata) {
    // find CELL_DATA, POINT_DATA which already present
    std::vector<size_t> fnd = find_strings_in_stream(fname, {"CELL_DATA", "POINT_DATA"});

    if (fnd[1] != std::string::npos) {
        // if point data exist write data before it
        std::ifstream f1(fname);
        std::string out(std::istreambuf_iterator<char>{f1}, std::istreambuf_iterator<char>{});
        f1.close();
        std::ofstream f2(fname);
        f2 << out.substr(0, fnd[1]);
        if (fnd[0] == std::string::npos)
            append_cell_data_header(ndata, f2);
        f2 << data << out.substr(fnd[1], out.npos);
    } else {
        // if not, simply append
        std::ofstream f2(fname, std::ios::app);
        if (fnd[0] == std::string::npos)
            append_cell_data_header(ndata, f2);
        f2 << data;
    }
}

void write_point_string(std::string fname, const std::string& data, size_t ndata) {
    // find if there are any point data already in file
    std::vector<size_t> fnd = find_strings_in_stream(fname, {"POINT_DATA"});

    // vertex data do not exist write caption
    std::fstream f2(fname, std::ios::app);
    if (fnd[0] == std::string::npos)
        append_point_data_header(ndata, f2);
    // write data
    f2 << data;
    f2.close();
}

} // namespace

void VtkUtils::append_header(std::string caption, std::ostream& fs) {
    fs << "# vtk DataFile Version 3.0" << std::endl;
    fs << caption << std::endl;
    fs << "ASCII" << std::endl;
}

void VtkUtils::append_points(const std::vector<Point>& points, std::ostream& fs) {
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
    fs << "POINTS " << points.size() << " double" << std::endl;
    for (const auto& point : points) {
        fs << point.x << " " << point.y << " " << point.z;
        fs << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////
// add cell info
///////////////////////////////////////////////////////////////////////////////

void VtkUtils::add_cell_data(const std::vector<double>& data, std::string data_cap, std::string fname, size_t ndata) {
    set_ndata(data, ndata);
    std::ostringstream fs;
    add_data_stream(data.begin(), ndata, data_cap, fs);
    write_cell_string(fname, fs.str(), ndata);
}

void VtkUtils::add_cell_vector(const std::vector<Vector>& data, std::string data_cap, std::string fname, size_t ndata) {
    set_ndata(data, ndata);
    std::ostringstream fs;
    add_vector_stream(data.begin(), ndata, data_cap, fs);
    write_cell_string(fname, fs.str(), ndata);
}

void VtkUtils::add_cell_vector(const std::vector<double>& ux, const std::vector<double>& uy, std::string data_cap,
                               std::string fname, size_t ndata) {
    set_ndata(ux, ndata);
    std::ostringstream fs;
    add_vector_stream(ux.begin(), uy.begin(), ndata, data_cap, fs);
    write_cell_string(fname, fs.str(), ndata);
}

void VtkUtils::add_cell_vector(const std::vector<double>& ux, const std::vector<double>& uy,
                               const std::vector<double>& uz, std::string data_cap, std::string fname, size_t ndata) {

    set_ndata(ux, ndata);
    std::ostringstream fs;
    add_vector_stream(ux.begin(), uy.begin(), uz.begin(), ndata, data_cap, fs);
    write_cell_string(fname, fs.str(), ndata);
}

///////////////////////////////////////////////////////////////////////////////
// add point info
///////////////////////////////////////////////////////////////////////////////

void VtkUtils::add_point_data(const std::vector<double>& data, std::string data_cap, std::string fname, size_t ndata) {
    set_ndata(data, ndata);
    std::ostringstream fs;
    add_data_stream(data.begin(), ndata, data_cap, fs);
    write_point_string(fname, fs.str(), ndata);
}

void VtkUtils::add_point_vector(const std::vector<Vector>& data, std::string data_cap, std::string fname,
                                size_t ndata) {
    set_ndata(data, ndata);
    std::ostringstream fs;
    add_vector_stream(data.begin(), ndata, data_cap, fs);
    write_point_string(fname, fs.str(), ndata);
}

void VtkUtils::add_point_vector(const std::vector<double>& ux, const std::vector<double>& uy, std::string data_cap,
                                std::string fname, size_t ndata) {
    set_ndata(ux, ndata);
    std::ostringstream fs;
    add_vector_stream(ux.begin(), uy.begin(), ndata, data_cap, fs);
    write_point_string(fname, fs.str(), ndata);
}

void VtkUtils::add_point_vector(const std::vector<double>& ux, const std::vector<double>& uy,
                                const std::vector<double>& uz, std::string data_cap, std::string fname, size_t ndata) {
    set_ndata(ux, ndata);
    std::ostringstream fs;
    add_vector_stream(ux.begin(), uy.begin(), uz.begin(), ndata, data_cap, fs);
    write_point_string(fname, fs.str(), ndata);
}

///////////////////////////////////////////////////////////////////////////////
// TimeSeriesWriter
///////////////////////////////////////////////////////////////////////////////

VtkUtils::TimeSeriesWriter::TimeSeriesWriter(std::string stem) : stem_(stem), series_fn_(stem + ".vtk.series") {
    create_directory(stem, true);

    std::ofstream ofs(series_fn_);
    if (!ofs)
        throw std::runtime_error("Failed to open " + series_fn_ + " for writing");

    ofs << "{" << std::endl;
    ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
    ofs << "  \"files\" : [" << std::endl;
    ofs << "  ]" << std::endl;
    ofs << "}" << std::endl;
}

std::string VtkUtils::TimeSeriesWriter::add(double tm) {
    if (step_ > 0) {
        int itime_point = get_time_point_index(tm);
        if (itime_point > last_saved_point_) {
            last_saved_point_ = itime_point;
        } else {
            return "";
        }
    }

    std::ostringstream fn;
    fn << std::setfill('0') << std::setw(8) << std::fixed << std::setprecision(4) << tm << ".vtk";
    std::string ret = stem_ + '/' + fn.str();

    std::ostringstream oss;
    if (!fileslist_.empty()) {
        oss << "," << std::endl;
    }
    oss << "    {\"name\": \"" << ret << "\", \"time\": " << tm << "}";
    fileslist_ += oss.str();

    save_series();

    return ret;
}

void VtkUtils::TimeSeriesWriter::save_series() const {
    std::fstream ofs(series_fn_);
    if (!ofs)
        throw std::runtime_error("Failed to open " + series_fn_ + " for writing");

    ofs << "{" << std::endl;
    ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
    ofs << "  \"files\" : [" << std::endl;
    ofs << fileslist_ << std::endl;
    ofs << "  ]" << std::endl;
    ofs << "}" << std::endl;
}

int VtkUtils::TimeSeriesWriter::get_time_point_index(double tm) const {
    int index = int(std::floor(tm / step_));
    double cur_point = index * step_;
    double next_point = cur_point + step_;
    if (next_point - tm < step_eps_) {
        return index + 1;
    } else {
        return index;
    }
}

void VtkUtils::TimeSeriesWriter::set_time_step(double tm_step, double eps) {
    step_ = tm_step;
    step_eps_ = eps;
}
