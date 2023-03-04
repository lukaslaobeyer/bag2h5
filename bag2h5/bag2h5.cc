#include "logging.h"

#include <vector>
#include <unordered_map>

#include <Eigen/Dense>
#include <magic_enum.hpp>
#include <embag/view.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>


namespace {

namespace h5 = HighFive;
using NDIndex = std::vector<size_t>;


template <typename T>
std::vector<T> append(const std::vector<T>& vec, const T& val) {
    std::vector<size_t> new_vec(vec.size() + 1);
    std::copy(vec.begin(), vec.end(), new_vec.begin());
    new_vec[vec.size()] = val;
    return new_vec;
}

template <typename T>
std::string it_to_string(const T& vals, const std::string& sep) {
    std::string str("");
    for (const auto& v : vals) {
        str.append(std::to_string(v) + sep);
    }
    str.pop_back();
    return str;
}

size_t get_size(const Embag::RosValue::Pointer& v) {
    if (v->getType() == Embag::RosValue::Type::primitive_array) {
        return v->getPrimitiveArrayRosValueBufferSize();
    }

    if (v->getType() == Embag::RosValue::Type::string) {
        return v->as<uint32_t>();
    }

    return Embag::RosValue::primitiveTypeToSize(v->getType());
}

bool is_primitive(const Embag::RosValue::Type& t) {
    switch (t) {
    case Embag::RosValue::Type::ros_bool:
    case Embag::RosValue::Type::int8:
    case Embag::RosValue::Type::uint8:
    case Embag::RosValue::Type::int16:
    case Embag::RosValue::Type::uint16:
    case Embag::RosValue::Type::int32:
    case Embag::RosValue::Type::uint32:
    case Embag::RosValue::Type::int64:
    case Embag::RosValue::Type::uint64:
    case Embag::RosValue::Type::float32:
    case Embag::RosValue::Type::float64:
    case Embag::RosValue::Type::string:
    case Embag::RosValue::Type::ros_time:
    case Embag::RosValue::Type::ros_duration:
    case Embag::RosValue::Type::primitive_array:
        return true;

    case Embag::RosValue::Type::object:
    case Embag::RosValue::Type::array:
        return false;

    default:
        throw std::runtime_error("should never happen");
    }
}

h5::DataType h5_datatype_from_embag_type(const Embag::RosValue::Type& t) {
    switch (t) {
    case Embag::RosValue::Type::ros_bool:
        return h5::create_datatype<bool>();
    case Embag::RosValue::Type::int8:
        return h5::create_datatype<int8_t>();
    case Embag::RosValue::Type::uint8:
        return h5::create_datatype<uint8_t>();
    case Embag::RosValue::Type::int16:
        return h5::create_datatype<int16_t>();
    case Embag::RosValue::Type::uint16:
        return h5::create_datatype<uint16_t>();
    case Embag::RosValue::Type::int32:
        return h5::create_datatype<int32_t>();
    case Embag::RosValue::Type::uint32:
        return h5::create_datatype<uint32_t>();
    case Embag::RosValue::Type::int64:
        return h5::create_datatype<int64_t>();
    case Embag::RosValue::Type::uint64:
        return h5::create_datatype<uint64_t>();
    case Embag::RosValue::Type::float32:
        return h5::create_datatype<float>();
    case Embag::RosValue::Type::float64:
        return h5::create_datatype<double>();
    case Embag::RosValue::Type::string:
        return h5::create_datatype<std::string>();
    case Embag::RosValue::Type::ros_time:
        return h5::create_datatype<double>();
    case Embag::RosValue::Type::ros_duration:
        return h5::create_datatype<double>();

    case Embag::RosValue::Type::primitive_array:
    case Embag::RosValue::Type::object:
    case Embag::RosValue::Type::array:
        throw std::runtime_error("only primitive and non-array type supported");

    default:
        throw std::runtime_error("should never happen");
    }
}

h5::DataType h5_datatype_from_ros_value(const Embag::RosValue::Pointer& v) {
    const Embag::RosValue::Type t = v->getType();
    if (t == Embag::RosValue::Type::primitive_array)
        return h5_datatype_from_embag_type(v->getElementType());
    return h5_datatype_from_embag_type(t);
}

std::vector<hsize_t> get_primitive_array_dims(const Embag::RosValue::Pointer& v) {
    const Embag::RosValue::Type t = v->getType();
    if (!is_primitive(t))
        throw std::runtime_error("input value must be primitive type");
    if (t == Embag::RosValue::Type::primitive_array)
        return {v->size()};
    return {};
}

struct H5Dataset {
    H5Dataset(h5::DataSet&& d_in,
              std::vector<hsize_t>&& chunk_dims_in,
              std::vector<hsize_t>&& prim_dims_in,
              const Embag::RosValue::Type& type_in,
              const h5::DataType& h5_type_in)
        : d(d_in), chunk_dims(chunk_dims_in), prim_dims(prim_dims_in), type(type_in), h5_type(h5_type_in) {}

    h5::DataSet d;
    const std::vector<hsize_t> chunk_dims;
    const std::vector<hsize_t> prim_dims;
    const Embag::RosValue::Type type;
    const h5::DataType h5_type;
};

H5Dataset create_dataset(h5::File& f,
                           const std::string& dataset_name,
                           const NDIndex& dims,
                           const Embag::RosValue::Pointer& v) {
    // Dataspace
    std::vector<hsize_t> primitive_array_dims = get_primitive_array_dims(v);
    std::vector<hsize_t> full_dims(1 + dims.size() + primitive_array_dims.size());
    full_dims.at(0) = 0;
    std::transform(dims.begin(),
                   dims.end(),
                   full_dims.begin() + 1,
                   [](const size_t& i) { return static_cast<hsize_t>(i); });
    std::copy(primitive_array_dims.begin(), primitive_array_dims.end(), full_dims.begin() + 1 + dims.size());
    std::vector<hsize_t> max_dims(full_dims);
    max_dims[0] = h5::DataSpace::UNLIMITED;
    h5::DataSpace dataspace(full_dims, max_dims);

    // Chunking
    constexpr hsize_t MAX_FLAT_CHUNK_SIZE = 1000000;

    const size_t nonprim_dims_prod = std::reduce(
        dims.begin(), dims.end(), 1, std::multiplies<>()
    );
    const size_t elem_size = nonprim_dims_prod * get_size(v);
    std::vector<hsize_t> chunk_dims(full_dims.size());
    chunk_dims.at(0) = 1 + std::trunc(static_cast<double>(MAX_FLAT_CHUNK_SIZE)/static_cast<double>(elem_size) - 0.5);
    std::copy(full_dims.begin() + 1, full_dims.end(), chunk_dims.begin() + 1);

    // Datatype
    h5::DataType datatype = h5_datatype_from_ros_value(v);

    // Debug
    LOG_DEBUG("Creating dataspace for %s [%s]", dataset_name.c_str(), it_to_string(dims, "x").c_str());
    LOG_DEBUG("num prim dims %lu", primitive_array_dims.size());
    LOG_DEBUG("dataspace: %s, chunks: %s", it_to_string(full_dims, "x").c_str(), it_to_string(chunk_dims, "x").c_str());

    h5::DataSetCreateProps props;
    props.add(h5::Chunking(chunk_dims));
    return H5Dataset(
        f.createDataSet(dataset_name, dataspace, datatype, props),
        std::move(chunk_dims),
        std::move(primitive_array_dims),
        v->getType(),
        datatype
    );
}

class NDBuffer {
public:
    NDBuffer(h5::File& f,
             const std::string& dataset_name,
             const NDIndex& dims,
             const Embag::RosValue::Pointer& v)
        : d_(create_dataset(f, dataset_name, dims, v))
        , dims_(dims)
        , awaiting_first_msg_(true)
        , dataset_index_(d_.chunk_dims.size(), 0)
        , dataset_index_count_(d_.chunk_dims.size(), 1)
        , dataset_size_(d_.chunk_dims)
        , msg_idx_(0) {
        std::copy(d_.prim_dims.begin(),
                  d_.prim_dims.end(),
                  dataset_index_count_.end() - d_.prim_dims.size());
        if (d_.prim_dims.size() == 0 && d_.type != Embag::RosValue::Type::string) {
            // use the buffer except for primitive array type
            buf_.resize(get_size(v));
        }
    }

    const NDIndex& dims() const {
        return dims_;
    }

    void insert(const size_t msg_idx, const NDIndex& index, const Embag::RosValue::Pointer& v) {
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

        // Copy into dataset
        std::copy(index.begin(), index.end(), dataset_index_.begin() + 1);
        if (d_.type == Embag::RosValue::Type::primitive_array) {
            d_.d.select(dataset_index_, dataset_index_count_).write_raw(
                reinterpret_cast<const char*>(v->getPrimitiveArrayRosValueBuffer()),
                d_.h5_type
            );
        } else if (d_.type == Embag::RosValue::Type::string) {
            const std::string str = v->as<std::string>();
            d_.d.select(dataset_index_, dataset_index_count_).write(str);
        } else if (d_.type == Embag::RosValue::Type::ros_time) {
            const double timestamp = v->as<Embag::RosValue::ros_time_t>().to_sec();
            d_.d.select(dataset_index_, dataset_index_count_).write_raw(
                &timestamp, d_.h5_type
            );
        } else if (d_.type == Embag::RosValue::Type::ros_duration) {
            const double duration = v->as<Embag::RosValue::ros_duration_t>().to_sec();
            d_.d.select(dataset_index_, dataset_index_count_).write_raw(
                &duration, d_.h5_type
            );
        } else {
            d_.d.select(dataset_index_, dataset_index_count_).write_raw(
                v->data(), d_.h5_type
            );
        }
    }

private:
    H5Dataset d_;
    const NDIndex dims_;

    std::vector<char> buf_;

    bool awaiting_first_msg_;
    std::vector<hsize_t> dataset_index_;
    std::vector<hsize_t> dataset_index_count_;
    std::vector<hsize_t> dataset_size_;
    size_t msg_idx_;
};

class H5Writer {
public:
    H5Writer(std::string outfile)
        : h5_file_(outfile, h5::File::Overwrite) {}

    void process_msg(Embag::RosMessage& m) {
        process_value(num_msgs_++, m.data(), m.topic);
    }

    void print_datasets() const {
        for (const auto& [key, val] : datasets_) {
            if (val.dims().size() > 0)
                LOG_INFO("%s [%s]", key.c_str(), it_to_string(val.dims(), "x").c_str());
            else
                LOG_INFO("%s", key.c_str());
        }
    }

private:
    void process_value(const size_t msg_idx,
                       const Embag::RosValue::Pointer& v,
                       const std::string& dataset="",
                       const NDIndex& dims={},
                       const NDIndex& index={}) {
        if (is_primitive(v->getType())) {
            const std::string type_name(magic_enum::enum_name(v->getType()));
            insert_value(msg_idx, dataset, dims, index, v);
        } else if (v->getType() == Embag::RosValue::Type::object) {
            for (const auto& [key, child] : v->getObjects()) {
                process_value(msg_idx, child, dataset + "/" + key, dims, index);
            }
        } else if (v->getType() == Embag::RosValue::Type::array) {
            const auto& values = v->getValues();
            if (values.size() > 0) {
                NDIndex new_dims(append(dims, values.size()));
                for (size_t i = 0; i < values.size(); ++i) {
                    NDIndex new_index(append(index, i));
                    process_value(msg_idx,
                                  values.at(i),
                                  dataset,
                                  new_dims,
                                  new_index);
                }
            }
        }
    }

    void insert_value(const size_t msg_idx,
                      const std::string& dataset,
                      const NDIndex& dims,
                      const NDIndex& index,
                      const Embag::RosValue::Pointer& v) {
        if (!is_primitive(v->getType())) {
            throw std::runtime_error("cannot write non-primitive type");
        }

        if (datasets_.count(dataset) == 0) {
            datasets_.emplace(std::piecewise_construct,
                              std::forward_as_tuple(dataset),
                              std::forward_as_tuple(h5_file_, dataset, dims, v));
        }

        NDBuffer& d = datasets_.at(dataset);
        if (d.dims() != dims) {
            throw std::runtime_error("dims must be constant for each field");
        }
        d.insert(msg_idx, index, v);
    }

private:
    size_t num_msgs_;
    h5::File h5_file_;
    std::unordered_map<std::string, NDBuffer> datasets_;
};

}

int main(int argc, char** argv) {
    if (argc < 3) {
        LOG_ERROR("usage: file1.bag [file2.bag ...] out.h5");
        return -1;
    }

    std::vector<std::string> bagfiles(argv + 1, argv + argc - 1);
    std::string outfile(argv[argc - 1]);

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

    H5Writer h5w(outfile);

    for (const auto& m : view.getMessages()) {
        h5w.process_msg(*m);

        const double t = m->timestamp.to_sec() - t_start;
        if (std::isnan(t_prev) || t - t_prev > log_dt) {
            logging::print_progress(t/(t_end - t_start));
            t_prev = t;
        }
    }
    logging::print_progress(1.0);

    h5w.print_datasets();

    return 0;
}
