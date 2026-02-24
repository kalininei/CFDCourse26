#ifndef __CFD_GEOM_VTK_HPP__
#define __CFD_GEOM_VTK_HPP__

#include "cfd/grid/i_grid.hpp"
#include <fstream>

namespace cfd {

/**
 * @brief Collection of vtk-writer utilities
 */
struct VtkUtils {
    /// adds vtk caption
    static void append_header(std::string caption, std::ostream& fs);

    /// adds point list
    static void append_points(const std::vector<Point>& points, std::ostream& fs);

    /**
     * @brief adds cell data to saved vtk grid
     * @param data      data vector
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     file name of saved vtk grid
     * @param ndata     length of data vector. If =0 uses data.size() as its
     * length
     */
    static void add_cell_data(const std::vector<double>& data, std::string data_cap, std::string fname,
                              size_t ndata = 0);

    /**
     * @brief adds cell vector data to saved vtk file
     * @param data      data vector
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses data.size() as its
     * length
     */
    static void add_cell_vector(const std::vector<Vector>& data, std::string data_cap, std::string fname,
                                size_t ndata = 0);

    /**
     * @brief adds 2d cell vector data to saved vtk file
     * @param ux        vector x coordinates
     * @param uy        vector y coordinates
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses ux.size() as its length
     */
    static void add_cell_vector(const std::vector<double>& ux, const std::vector<double>& uy, std::string data_cap,
                                std::string fname, size_t ndata = 0);

    /**
     * @brief adds 3d cell vector data to saved vtk file
     * @param ux        vector x coordinates
     * @param uy        vector y coordinates
     * @param uz        vector z coordinates
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses ux.size() as its length
     */
    static void add_cell_vector(const std::vector<double>& ux, const std::vector<double>& uy,
                                const std::vector<double>& uz, std::string data_cap, std::string fname,
                                size_t ndata = 0);

    /**
     * @brief adds vertex data to the saved vtk file
     * @param data      data vector
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses data.size() as its
     * length
     */
    static void add_point_data(const std::vector<double>& data, std::string data_cap, std::string fname,
                               size_t ndata = 0);

    /**
     * @brief adds vertex vector data to saved vtk file
     * @param data      data vector
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses data.size() as its
     * length
     */
    static void add_point_vector(const std::vector<Vector>& data, std::string data_cap, std::string fname,
                                 size_t ndata = 0);

    /**
     * @brief adds 2d point vector data to saved vtk file
     * @param ux        vector x coordinates
     * @param uy        vector y coordinates
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses ux.size() as its length
     */
    static void add_point_vector(const std::vector<double>& ux, const std::vector<double>& uy, std::string data_cap,
                                 std::string fname, size_t ndata = 0);

    /**
     * @brief adds 3d point vector data to saved vtk file
     * @param ux        vector x coordinates
     * @param uy        vector y coordinates
     * @param uz        vector z coordinates
     * @param data_cap  data caption as it should appear in vtk
     * @param fname     saved vtk file name
     * @param ndata     length of data vector. If =0 uses ux.size() as its length
     */
    static void add_point_vector(const std::vector<double>& ux, const std::vector<double>& uy,
                                 const std::vector<double>& uz, std::string data_cap, std::string fname,
                                 size_t ndata = 0);

    /**
     * @brief Writes data for increasing time points into the .vtk.series file
     */
    class TimeSeriesWriter {
    public:
        /**
         * @brief constructor
         * @param stem file name stem
         *
         * This creates the series file named "<stem>.vtk.series".
         * Vtk files with instant data will be saved into the "<stem>" directory.
         * If this directory already exists it will be purged.
         */
        TimeSeriesWriter(std::string stem);

        /**
         * @brief adds new time point to the series file
         * @param   tm time point value
         * @returns file name in <stem> directory or empty string if the given time
         * point is not valid for save due to the time step condition
         *
         * This function only makes a record in the series file but does not create
         * vtk file with instant data in the "<stem>" directory. The latter should
         * be done manually using the returned filename.
         */
        std::string add(double tm, bool force = false);
        std::string add_iter(size_t iter, bool force = false);

        /**
         * @breif set saving time step
         * @param tm_step step value
         *
         * Default step value is zero, that will make saving for each time point.
         * If non zero value is set as a step then time save points will be
         * calculated only a single time point will be chosen for [N*tm_step-eps,
         * (N+1)*tm_step] time period. TimeDependentWriter::add() calls with ignored
         * time points will return empty strings.
         */
        void set_time_step(double tm_step, double eps = 1e-6);
        void set_iter_step(size_t iter_step);

    private:
        const std::string stem_;
        const std::string series_fn_;
        std::string fileslist_;
        double step_ = 0;
        double step_eps_ = 0;
        int last_saved_point_ = -1;
        void save_series() const;
        int get_time_point_index(double tm) const;
    };
};

} // namespace cfd

#endif
