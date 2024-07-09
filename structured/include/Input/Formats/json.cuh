#pragma once

#define JSON_TYPE nlohmann::json

inline void to_json(JSON_TYPE& j, const float2& p) {
    j = JSON_TYPE{p.x,p.y};
}

inline void from_json(const JSON_TYPE& j, float2& p) {
    p.x=float(j[0]);
    p.y=float(j[1]);
}


inline void to_json(JSON_TYPE& j, const float3& p) {
    j = JSON_TYPE{p.x,p.y,p.z};
}

inline void from_json(const JSON_TYPE& j, float3& p) {
    p.x=float(j[0]);
    p.y=float(j[1]);
    p.z=float(j[2]);
}

inline void to_json(JSON_TYPE& j, const float4& p) {
    j = JSON_TYPE{p.x,p.y,p.z,p.w};
}

inline void from_json(const JSON_TYPE& j, float4& p) {
    p.x=float(j[0]);
    p.y=float(j[1]);
    p.z=float(j[2]);
    p.w=float(j[3]);
}

inline void to_json(JSON_TYPE& j, const double2& p) {
    j = JSON_TYPE{p.x,p.y};
}

inline void from_json(const JSON_TYPE& j, double2& p) {
    p.x=double(j[0]);
    p.y=double(j[1]);
}

inline void to_json(JSON_TYPE& j, const double3& p) {
    j = JSON_TYPE{p.x,p.y,p.z};
}

inline void from_json(const JSON_TYPE& j, double3& p) {
    p.x=double(j[0]);
    p.y=double(j[1]);
    p.z=double(j[2]);
}

inline void to_json(JSON_TYPE& j, const double4& p) {
    j = JSON_TYPE{p.x,p.y,p.z,p.w};
}

inline void from_json(const JSON_TYPE& j, double4& p) {
    p.x=double(j[0]);
    p.y=double(j[1]);
    p.z=double(j[2]);
    p.w=double(j[3]);
}

inline nlohmann::json_pointer<nlohmann::json::basic_json::string_t>
path2jsonPointer(const std::vector<std::string>& path){
    std::string path_ptr = "";

    for(const std::string& p : path){
        path_ptr+="/"+p;
    }

    return nlohmann::json_pointer<std::string>(path_ptr);
}

