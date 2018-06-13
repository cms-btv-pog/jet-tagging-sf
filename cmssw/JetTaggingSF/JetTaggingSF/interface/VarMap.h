/*
 * Objects to handle and store variables.
 *
 * Authors:
 *   - Marcel Rieger
 */

#ifndef JETTAGGINGSF_JETTAGGINGSF_H
#define JETTAGGINGSF_JETTAGGINGSF_H

#include <algorithm>

typedef std::string string;
typedef std::vector<string> vstring;

class VarMap
{
public:
    VarMap()
    {
    }

    ~VarMap()
    {
    }

    inline size_t size() const
    {
        return names_.size();
    }

    inline string getName(size_t i) const
    {
        return names_[i];
    }

    inline bool contains(const string name) const
    {
        return std::find(names_.begin(), names_.end(), name) != names_.end();
    }

    inline string getFlag(const string name)
    {
        return typeFlags_[name];
    }

    void reset()
    {
        std::map<string, float>::iterator itFloat;
        for (itFloat = floatValues_.begin(); itFloat != floatValues_.end(); itFloat++)
        {
            itFloat->second = emptyFloat;
        }

        std::map<string, double>::iterator itDouble;
        for (itDouble = doubleValues_.begin(); itDouble != doubleValues_.end(); itDouble++)
        {
            itDouble->second = emptyDouble;
        }

        std::map<string, int32_t>::iterator itInt32;
        for (itInt32 = int32Values_.begin(); itInt32 != int32Values_.end(); itInt32++)
        {
            itInt32->second = emptyInt32;
        }

        std::map<string, int64_t>::iterator itInt64;
        for (itInt64 = int64Values_.begin(); itInt64 != int64Values_.end(); itInt64++)
        {
            itInt64->second = emptyInt64;
        }
    }

    void addFloat(const string name)
    {
        complainOnDuplicate(name);
        names_.push_back(name);
        floatValues_[name] = emptyFloat;
        typeFlags_[name] = "F";
    }

    inline void setFloat(const string name, float value)
    {
        floatValues_[name] = value;
    }

    inline float& getFloat(const string name)
    {
        return floatValues_[name];
    }

    void addDouble(const string name)
    {
        complainOnDuplicate(name);
        names_.push_back(name);
        doubleValues_[name] = emptyDouble;
        typeFlags_[name] = "D";
    }

    inline void setDouble(const string name, double value)
    {
        doubleValues_[name] = value;
    }

    inline double& getDouble(const string name)
    {
        return doubleValues_[name];
    }

    void addInt32(const string name)
    {
        complainOnDuplicate(name);
        names_.push_back(name);
        int32Values_[name] = emptyInt32;
        typeFlags_[name] = "I";
    }

    inline void setInt32(const string name, int32_t value)
    {
        int32Values_[name] = value;
    }

    inline int32_t& getInt32(const string name)
    {
        return int32Values_[name];
    }

    void addInt64(const string name)
    {
        complainOnDuplicate(name);
        names_.push_back(name);
        int64Values_[name] = emptyInt64;
        typeFlags_[name] = "L";
    }

    inline void setInt64(const string name, int64_t value)
    {
        int64Values_[name] = value;
    }

    inline int64_t& getInt64(const string name)
    {
        return int64Values_[name];
    }

private:
    float emptyFloat = -1e5;
    double emptyDouble = -1e5;
    int32_t emptyInt32 = int32_t(-1e5);
    int64_t emptyInt64 = int64_t(-1e5);

    vstring names_;
    std::map<string, string> typeFlags_;
    std::map<string, float> floatValues_;
    std::map<string, double> doubleValues_;
    std::map<string, int32_t> int32Values_;
    std::map<string, int64_t> int64Values_;

    inline void complainOnDuplicate(const string name) const
    {
        if (contains(name))
        {
            throw std::runtime_error("duplicate variable " + name);
        }
    }
};

#endif // JETTAGGINGSF_JETTAGGINGSF_H
