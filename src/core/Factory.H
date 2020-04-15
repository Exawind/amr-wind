#ifndef FACTORY_H
#define FACTORY_H

#include <memory>
#include <unordered_map>
#include <iostream>

namespace amr_wind {

template<class Base, class ... Args>
struct Factory
{
    static std::unique_ptr<Base>
    create(const std::string& key, Args... args)
    {
        // TODO: Add error diagnostics
        return table().at(key)(std::forward<Args>(args)...);
    }

    static void print(std::ostream& os)
    {
        const auto& tbl = table();
        os << Base::base_identifier() << " " << tbl.size() << std::endl;
        for (const auto& it: tbl) {
            os << " - " << it.first << std::endl;
        }
    }

    template<class T>
    struct Register : public Base
    {
        friend T;
        static bool registered;

        static bool add_sub_type()
        {
            // TODO: Add checks against multiple registration
            Factory::table()[T::identifier()] =
                [](Args... args) -> std::unique_ptr<Base> {
                return std::unique_ptr<Base>(
                    new T(std::forward<Args>(args)...));
            };
            return true;
        }

        virtual ~Register()
        {
            if (registered) {
                auto& tbl = Factory::table();
                const auto& it = tbl.find(T::identifier());
                if (it != tbl.end()) {
                    tbl.erase(T::identifier());
                    registered = false;
                }
            }
        }

    private:
        Register() { (void)registered; }
    };

    virtual ~Factory() = default;
    friend Base;

private:
    using CreatorFunc = std::unique_ptr<Base> (*)(Args...);
    using LookupTable = std::unordered_map<std::string, CreatorFunc>;

    static LookupTable& table()
    {
        static LookupTable tbl;
        return tbl;
    }

    Factory() = default;
};

template<class Base, class... Args>
template<class T>
bool Factory<Base, Args...>::Register<T>::registered =
    Factory<Base, Args...>::Register<T>::add_sub_type();

}

#endif /* FACTORY_H */
