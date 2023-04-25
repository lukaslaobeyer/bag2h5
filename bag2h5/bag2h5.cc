#include "logging.h"
#include "type_info.h"

#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <cxxopts.hpp>
#include <embag/view.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <magic_enum.hpp>

namespace bag2h5::h5_conversion {

namespace h5 = HighFive;
using NDIndex = std::vector<size_t>;

template <typename T>
std::vector<T> append(const std::vector<T>& vec, const T& val) {
    std::vector<T> new_vec(vec.size() + 1);
    std::copy(vec.begin(), vec.end(), new_vec.begin());
    new_vec[vec.size()] = val;
    return new_vec;
}

template <typename T>
std::vector<T> append(const std::vector<T>& vec0, const std::vector<T>& vec1) {
    std::vector<T> new_vec(vec0.size() + vec1.size());
    std::copy(vec0.begin(), vec0.end(), new_vec.begin());
    std::copy(vec1.begin(), vec1.end(), new_vec.begin() + vec0.size());
    return new_vec;
}

template <typename T>
std::string it_to_string(const T& vals, const std::string& sep) {
    if (vals.size() == 0) return "<empty>";

    std::string str("");
    for (const auto& v : vals) {
        str.append(std::to_string(v) + sep);
    }
    str.pop_back();
    return str;
}

class DatasetInfo {
public:
    DatasetInfo(const NDIndex& dims, const NDIndex& prim_dims, const h5::DataType& h5_type_in)
        : h5_type(h5_type_in)
        , dims_(dims)
        , prim_dims_(prim_dims)
        , is_variable_length_(false) {}

    void resize(const NDIndex& new_dims, const NDIndex& new_prim_dims) {
        if (new_dims.size() != dims_.size() || new_prim_dims.size() != prim_dims_.size())
            throw new std::runtime_error("variable dimensionality not supported");

        for (size_t i = 0; i < dims_.size(); i++) {
            if (dims_.at(i) != new_dims.at(i)) {
                is_variable_length_ = true;
                if (new_dims.at(i) > dims_.at(i)) dims_.at(i) = new_dims.at(i);
            }
        }

        for (size_t i = 0; i < prim_dims_.size(); i++) {
            if (prim_dims_.at(i) != new_prim_dims.at(i)) {
                is_variable_length_ = true;
                if (new_prim_dims.at(i) > prim_dims_.at(i)) prim_dims_.at(i) = new_prim_dims.at(i);
            }
        }
    }

    const NDIndex& dims() const { return dims_; }
    const NDIndex& prim_dims() const { return prim_dims_; }
    bool is_variable_length() const { return is_variable_length_; }

public:
    const h5::DataType h5_type;

private:
    NDIndex dims_;
    NDIndex prim_dims_;
    bool is_variable_length_;
};

struct H5Dataset {
    H5Dataset(h5::DataSet&& d_in,
              std::vector<hsize_t>&& chunk_dims_in,
              const std::vector<hsize_t>& prim_dims_in,
              const h5::DataType& h5_type_in)
        : d(d_in)
        , chunk_dims(chunk_dims_in)
        , prim_dims(prim_dims_in)
        , h5_type(h5_type_in) {}

    h5::DataSet d;
    const std::vector<hsize_t> chunk_dims;
    const std::vector<hsize_t> prim_dims;
    const h5::DataType h5_type;
};

class ROSDatasetAttributes {
public:
    ROSDatasetAttributes(Embag::RosMessage& m)
        : topic_name(m.topic)
        , msg_type(m.getTypeName()) {}

    void write(h5::DataSet& dataset) const {
        dataset.createAttribute<std::string>("ros_topic", h5::DataSpace::From(topic_name))
            .write(topic_name);
        dataset.createAttribute<std::string>("ros_msg_type", h5::DataSpace::From(msg_type))
            .write(msg_type);
    }

private:
    const std::string topic_name;
    const std::string msg_type;
};

H5Dataset create_dataset(h5::File& f,
                         const ROSDatasetAttributes& attrs,
                         const std::string& dataset_name,
                         const DatasetInfo& info) {
    using namespace type_info;

    // Dataspace
    std::vector<hsize_t> full_dims(1 + info.dims().size() + info.prim_dims().size());
    full_dims.at(0) = 0;
    std::transform(info.dims().begin(), info.dims().end(), full_dims.begin() + 1,
                   [](const size_t& i) { return static_cast<hsize_t>(i); });
    std::copy(info.prim_dims().begin(), info.prim_dims().end(),
              full_dims.begin() + 1 + info.dims().size());
    std::vector<hsize_t> max_dims(full_dims);
    max_dims[0] = h5::DataSpace::UNLIMITED;
    h5::DataSpace dataspace(full_dims, max_dims);

    // Chunking
    constexpr hsize_t MAX_FLAT_CHUNK_SIZE = 1000000;

    const size_t nonprim_dims_prod =
        std::reduce(info.dims().begin(), info.dims().end(), 1, std::multiplies<>());
    const size_t prim_dims_prod =
        std::reduce(info.prim_dims().begin(), info.prim_dims().end(), 1, std::multiplies<>());
    if (prim_dims_prod == 0)
        throw std::runtime_error("cannot create dataset with zero element size");
    const size_t elem_size =
        nonprim_dims_prod *
        std::max(1uz, prim_dims_prod * info.h5_type.getSize()); // must be positive
    std::vector<hsize_t> chunk_dims(full_dims.size());
    chunk_dims.at(0) =
        1 +
        std::trunc(static_cast<double>(MAX_FLAT_CHUNK_SIZE) / static_cast<double>(elem_size) - 0.5);
    std::copy(full_dims.begin() + 1, full_dims.end(), chunk_dims.begin() + 1);

    // Debug
    LOG_DEBUG("Creating dataspace for %s [%s]", dataset_name.c_str(),
              it_to_string(append(info.dims(), info.prim_dims()), "x").c_str());
    LOG_DEBUG("  - num primitive dims = %lu", info.prim_dims().size());
    LOG_DEBUG("  - dataspace: %s, chunks: %s", it_to_string(full_dims, "x").c_str(),
              it_to_string(chunk_dims, "x").c_str());

    h5::DataSetCreateProps props;
    props.add(h5::Chunking(chunk_dims));

    h5::DataSet h5_dataset = f.createDataSet(dataset_name, dataspace, info.h5_type, props);
    attrs.write(h5_dataset);
    return H5Dataset(std::move(h5_dataset), std::move(chunk_dims), info.prim_dims(), info.h5_type);
}

class NDBuffer {
public:
    NDBuffer(h5::File& f,
             const ROSDatasetAttributes& attrs,
             const std::string& dataset_name,
             const DatasetInfo& info)
        : d_(create_dataset(f, attrs, dataset_name, info))
        , awaiting_first_msg_(true)
        , dataset_index_(d_.chunk_dims.size(), 0)
        , dataset_index_count_(d_.chunk_dims.size(), 1)
        , dataset_size_(d_.chunk_dims)
        , msg_idx_(0) {}

    template <typename T>
    void write_value(const T& v, h5::Selection&& selection) {
        selection.write(v);
    }

    void write_value(const Embag::RosValue::ros_time_t& v, h5::Selection&& selection) {
        const double timestamp = v.to_sec();
        selection.write_raw(&timestamp, d_.h5_type);
    }

    void write_value(const Embag::RosValue::Pointer& v, h5::Selection&& selection) {
        const Embag::RosValue::Type type = v->getType();
        if (!type_info::is_primitive(type)) {
            throw std::runtime_error("cannot write non-primitive type");
        }

        if (type == Embag::RosValue::Type::primitive_array) {
            selection.write_raw(reinterpret_cast<const char*>(v->getPrimitiveArrayRosValueBuffer()),
                                d_.h5_type);
        } else if (type == Embag::RosValue::Type::string) {
            const std::string str = v->as<std::string>();
            selection.write(str);
        } else if (type == Embag::RosValue::Type::ros_time) {
            const double time = v->as<Embag::RosValue::ros_duration_t>().to_sec();
            selection.write_raw(&time, d_.h5_type);
        } else if (type == Embag::RosValue::Type::ros_duration) {
            const double duration = v->as<Embag::RosValue::ros_duration_t>().to_sec();
            selection.write_raw(&duration, d_.h5_type);
        } else {
            selection.write_raw(v->data(), d_.h5_type);
        }
    }

    template <typename ValueT>
    void insert(const size_t msg_idx, const NDIndex& index, const ValueT& v) {
        if (awaiting_first_msg_) {
            awaiting_first_msg_ = false;
            msg_idx_ = msg_idx;
        }
        if (msg_idx < msg_idx_) {
            throw std::runtime_error("msg_idx must be monotonically increasing");
        }
        if (msg_idx > msg_idx_) {
            msg_idx_ = msg_idx;
            dataset_index_.at(0) += 1;
        }

        // Resize dataset
        dataset_size_.at(0) = dataset_index_.at(0) + 1;
        d_.d.resize(dataset_size_);

        // Use correct primitive array dims
        if (d_.prim_dims.size() > 0) {
            const NDIndex prim_dims = type_info::get_primitive_array_dims(v);
            std::copy(prim_dims.begin(), prim_dims.end(),
                      dataset_index_count_.end() - prim_dims.size());
        }

        // Copy into dataset
        std::copy(index.begin(), index.end(), dataset_index_.begin() + 1);
        write_value(v, d_.d.select(dataset_index_, dataset_index_count_));
    }

    const std::vector<hsize_t>& size() const { return dataset_size_; }

private:
    H5Dataset d_;

    bool awaiting_first_msg_;
    std::vector<hsize_t> dataset_index_;
    std::vector<hsize_t> dataset_index_count_;
    std::vector<hsize_t> dataset_size_;
    size_t msg_idx_;
};

class H5Writer {
public:
    enum class Mode {
        SCAN,
        WRITE,
    };

public:
    H5Writer(const std::string& outfile)
        : mode_(Mode::SCAN)
        , h5_file_(outfile, h5::File::Overwrite) {}

    void scan_msg(Embag::RosMessage& m) {
        if (mode_ != Mode::SCAN) throw new std::runtime_error("cannot scan message: this is a bug");
        process_value(m, num_msgs_, m.data(), m.topic);
    }

    void process_msg(Embag::RosMessage& m) {
        if (mode_ != Mode::WRITE) mode_ = Mode::WRITE;
        handle_timestamp(m, num_msgs_, m.topic + "/__ros_msg/timestamp", {}, {}, m.timestamp);
        process_value(m, num_msgs_, m.data(), m.topic);
        ++num_msgs_;
    }

    void print_datasets() const {
        for (const auto& [key, d] : datasets_) {
            LOG_INFO("%s [%s]", key.c_str(), it_to_string(d.size(), "x").c_str());
        }
    }

    void print_scan_results() const {
        for (auto&& [name, info] : dataset_info_) {
            if (info.is_variable_length()) LOG_INFO("Variable length dataset: %s", name.c_str());
        }
    }

private:
    void process_value(Embag::RosMessage& m,
                       const size_t msg_idx,
                       const Embag::RosValue::Pointer& v,
                       const std::string& dataset = "",
                       const NDIndex& dims = {},
                       const NDIndex& index = {}) {
        if (v->getType() == Embag::RosValue::Type::object) {
            for (const auto& [key, child] : v->getObjects()) {
                process_value(m, msg_idx, child, dataset + "/" + key, dims, index);
            }
        } else if (v->getType() == Embag::RosValue::Type::array) {
            const auto& values = v->getValues();
            if (values.size() > 0) {
                NDIndex new_dims(append(dims, values.size()));
                for (size_t i = 0; i < values.size(); ++i) {
                    NDIndex new_index(append(index, i));
                    process_value(m, msg_idx, values.at(i), dataset, new_dims, new_index);
                }
            }
        } else if (v->getType() == Embag::RosValue::Type::ros_time) {
            handle_timestamp(m, msg_idx, dataset, dims, index,
                             v->as<Embag::RosValue::ros_time_t>());
        } else if (type_info::is_primitive(v->getType())) {
            const std::string type_name(magic_enum::enum_name(v->getType()));
            handle_value(m, msg_idx, dataset, dims, index, v);
        } else {
            throw std::runtime_error("unknown value at " + dataset);
        }
    }

    void handle_timestamp(Embag::RosMessage& m,
                          const size_t msg_idx,
                          const std::string& dataset,
                          const NDIndex& dims,
                          const NDIndex& index,
                          const Embag::RosValue::ros_time_t& stamp) {
        constexpr bool CONVERT_TIMESTAMP_TO_DOUBLE = false;
        if (CONVERT_TIMESTAMP_TO_DOUBLE) {
            handle_value(m, msg_idx, dataset, dims, index, stamp);
        } else {
            handle_value(m, msg_idx, dataset + "/secs", dims, index, stamp.secs);
            handle_value(m, msg_idx, dataset + "/nsecs", dims, index, stamp.nsecs);
        }
    }

    template <typename ValueT>
    void handle_value(Embag::RosMessage& m,
                      const size_t msg_idx,
                      const std::string& dataset,
                      const NDIndex& dims,
                      const NDIndex& index,
                      const ValueT& v) {
        switch (mode_) {
        case Mode::SCAN:
            return record_dataset_info(dataset, dims, v);
        case Mode::WRITE: {
            if (!dataset_info_.contains(dataset)) record_dataset_info(dataset, dims, v);
            return insert_value(m, msg_idx, dataset, dims, index, v);
        }
        }
    }

    template <typename ValueT>
    void insert_value(Embag::RosMessage& m,
                      const size_t msg_idx,
                      const std::string& dataset_in,
                      const NDIndex& dims,
                      const NDIndex& index,
                      const ValueT& v) {
        const DatasetInfo& info = dataset_info_.at(dataset_in);

        // If the dataset is variable length, then create a new group out of it
        // with two subkeys representing its data and size.
        const std::string dataset = dataset_in + (info.is_variable_length() ? "/data" : "");

        if (!datasets_.contains(dataset)) {
            // Skip fields whose content is empty.
            const size_t prim_dims_prod = std::reduce(
                info.prim_dims().begin(), info.prim_dims().end(), 1, std::multiplies<>());
            if (prim_dims_prod == 0) {
                LOG_DEBUG("skipping message at %s since it is empty (primitive dimensions %s)",
                          dataset.c_str(), it_to_string(info.prim_dims(), "x").c_str());
                return;
            }

            const ROSDatasetAttributes attrs(m);
            datasets_.emplace(std::piecewise_construct, std::forward_as_tuple(dataset),
                              std::forward_as_tuple(h5_file_, attrs, dataset, info));
        }

        NDBuffer& d = datasets_.at(dataset);
        const NDIndex prim_dims = type_info::get_primitive_array_dims(v);
        if (info.is_variable_length()) {
            // Write variable sizes as new datasets.
            const NDIndex full_dims = append(dims, prim_dims);
            if (std::all_of(index.cbegin(), index.cend(), [](size_t i) { return i == 0; })) {
                if (full_dims.size() == 1)
                    handle_value(m, msg_idx, dataset_in + "/size", {}, {}, full_dims.at(0));
                else
                    for (size_t i = 0; i < full_dims.size(); i++)
                        handle_value(m, msg_idx, dataset_in + "/size", {full_dims.size()}, {i},
                                     full_dims.at(i));
            }
        } else {
            if (info.dims() != dims) {
                throw std::runtime_error(
                    "dims must be constant for each field (offending dataset: " + dataset +
                    " was of size " + it_to_string(info.dims(), "x") + " but is now " +
                    it_to_string(dims, "x") + ")");
            }
            if (info.prim_dims() != prim_dims) {
                throw std::runtime_error(
                    "primitive array dims must be constant for each field (offending dataset: " +
                    dataset + " was of size " + it_to_string(info.prim_dims(), "x") +
                    " but is now " + it_to_string(prim_dims, "x") + ")");
            }
        }
        d.insert(msg_idx, index, v);
    }

    template <typename ValueT>
    void record_dataset_info(const std::string& dataset, const NDIndex& dims, const ValueT& v) {
        using namespace type_info;
        if (dataset_info_.count(dataset) == 0) {
            dataset_info_.emplace(std::piecewise_construct, std::forward_as_tuple(dataset),
                                  std::forward_as_tuple(dims, get_primitive_array_dims(v),
                                                        h5_datatype_from_value(v)));
        }
        dataset_info_.at(dataset).resize(dims, get_primitive_array_dims(v));
    }

private:
    Mode mode_;
    size_t num_msgs_;
    h5::File h5_file_;

    // Only used during Mode::SCAN:
    std::unordered_map<std::string, DatasetInfo> dataset_info_;

    // Only used during Mode::WRITE:
    std::unordered_map<std::string, NDBuffer> datasets_;
};

} // namespace bag2h5::h5_conversion

int main(int argc, char** argv) {
    cxxopts::Options options(argv[0], "Convert ros bags to HDF5.");
    // clang-format off
    options
        .positional_help("file1.bag [file2.bag...] out.h5")
        .add_options()
            ("h,help", "print this help message")
            ("e,exclude", "topics to exclude", cxxopts::value<std::vector<std::string>>())
            ("positional", "", cxxopts::value<std::vector<std::string>>());
    // clang-format on
    options.parse_positional({"positional"});
    auto parsed_options = options.parse(argc, argv);

    if (parsed_options.count("help") || parsed_options.count("positional") == 0 ||
        parsed_options["positional"].as<std::vector<std::string>>().size() < 2) {
        std::cout << options.help() << std::endl;
        return -1;
    }

    const auto exclude_topics = ([&]() {
        if (parsed_options.count("exclude")) {
            const auto topics_vec = parsed_options["exclude"].as<std::vector<std::string>>();
            return std::set<std::string>(topics_vec.begin(), topics_vec.end());
        }
        return std::set<std::string>();
    })();
    for (auto&& s : exclude_topics) LOG_DEBUG("exclude: %s", s.c_str());

    const auto in_out_filenames = parsed_options["positional"].as<std::vector<std::string>>();
    const std::vector<std::string> bagfiles(in_out_filenames.begin(), in_out_filenames.end() - 1);
    const std::string outfile(in_out_filenames.back());

    Embag::View view;
    for (const auto& filename : bagfiles) {
        try {
            view.addBag(filename);
        } catch (const std::runtime_error& err) {
            LOG_ERROR("file %s not found", filename.c_str());
            return -1;
        }
    }

    const double t_start = view.getStartTime().to_sec();
    const double t_end = view.getEndTime().to_sec();
    const double log_dt = std::min(20., (t_end - t_start) / 100.);
    double t_prev = std::numeric_limits<double>::quiet_NaN();

    bag2h5::h5_conversion::H5Writer h5w(outfile);

    // First pass: scan through without writing anything
    //   to detect variable length datasets
    LOG_INFO("Scanning dataset...");
    for (const auto& m : view.getMessages()) {
        if (exclude_topics.contains(m->topic)) continue;
        h5w.scan_msg(*m);
    }
    h5w.print_scan_results();

    LOG_INFO("Writing hdf5");
    for (const auto& m : view.getMessages()) {
        if (exclude_topics.contains(m->topic)) continue;

        h5w.process_msg(*m);

        const double t = m->timestamp.to_sec() - t_start;
        if (std::isnan(t_prev) || t - t_prev > log_dt) {
            logging::print_progress(t / (t_end - t_start));
            t_prev = t;
        }
    }
    logging::print_progress(1.0);

    h5w.print_datasets();

    return 0;
}
