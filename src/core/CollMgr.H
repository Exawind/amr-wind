#ifndef COLLMGR_H
#define COLLMGR_H

#include <string>
#include <memory>
#include <unordered_map>
#include "AMReX_Vector.H"

namespace amr_wind {

template<typename Collection, typename Type>
class CollMgr
{
public:
    using TypePtr = std::unique_ptr<Type>;
    using TypeVector = amrex::Vector<TypePtr>;

    CollMgr() = default;

    ~CollMgr() = default;

    template<class ... Args>
    Type& create(const std::string& key, Args ... args)
    {
        m_obj_vec.emplace_back(
            Type::create(key, std::forward<Args>(args)...));
        m_obj_map[key] = m_obj_vec.size();

        return *m_obj_vec.back();
    }

    TypeVector& objects() { return m_obj_vec; }

    bool contains(const std::string& key) const
    {
        auto it = m_obj_map.find(key);
        return (it != m_obj_map.end());
    }

    Type& operator()(const std::string& key)
    {
        return *m_obj_vec[m_obj_map.at(key)];
    }

    template<typename T>
    T& operator()(const std::string& key)
    {
        return dynamic_cast<T&>(operator()(key));
    }

protected:
    TypeVector m_obj_vec;

    std::unordered_map<std::string, int> m_obj_map;
};

}

#endif /* COLLMGR_H */
