#pragma once

#include <embag/ros_value.h>
#include <highfive/H5File.hpp>


namespace bag2h5::type_info {

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

HighFive::DataType h5_datatype_from_embag_type(const Embag::RosValue::Type& t) {
    namespace h5 = HighFive;

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

HighFive::DataType h5_datatype_from_value(const Embag::RosValue::ros_time_t&) {
    return h5_datatype_from_embag_type(Embag::RosValue::Type::ros_time);
}

HighFive::DataType h5_datatype_from_value(const Embag::RosValue::Pointer& v) {
    const Embag::RosValue::Type t = v->getType();
    if (t == Embag::RosValue::Type::primitive_array)
        return h5_datatype_from_embag_type(v->getElementType());
    return h5_datatype_from_embag_type(t);
}

template <typename T>
HighFive::DataType h5_datatype_from_value(const T&) {
    return HighFive::create_datatype<T>();
}

template <typename T>
size_t get_size(const T&) {
    return sizeof(T);
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

template <typename T>
std::vector<hsize_t> get_primitive_array_dims(const T&) {
    return {};
}

std::vector<hsize_t> get_primitive_array_dims(const Embag::RosValue::Pointer& v) {
    const Embag::RosValue::Type t = v->getType();
    if (!is_primitive(t))
        throw std::runtime_error("input value must be primitive type");
    if (t == Embag::RosValue::Type::primitive_array)
        return {v->size()};
    return {};
}

}
