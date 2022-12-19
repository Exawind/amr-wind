#include <string>
#include <algorithm>

#include "amr-wind/utilities/DerivedQuantity.H"
#include "amr-wind/utilities/io_utils.H"

namespace amr_wind {
namespace {

inline std::string strip_spaces(const std::string& inp)
{
    std::string str(inp);
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

std::pair<std::string, std::vector<std::string>>
parse_derived_qty(const std::string& key)
{
    auto popen = key.find("(");
    // If this is not a function type field then return it
    if (popen == std::string::npos) {
        return {key, {}};
    }

    auto pclose = key.find(")");
    if (pclose == std::string::npos) {
        amrex::Abort(
            "Error encountered when parsing derived field name: " + key);
    }

    // Get the function name
    auto dname = strip_spaces(key.substr(0, popen));

    // Collect all arguments required for the derived quantity function
    auto fargs = strip_spaces(key.substr(popen + 1, pclose - popen - 1));
    if (fargs.back() != ',') {
        fargs += ",";
    }
    std::vector<std::string> args;
    size_t start = 0;
    size_t pos;
    while ((pos = fargs.find(",", start)) != std::string::npos) {
        args.push_back(fargs.substr(start, pos - start));
        start = pos + 1;
    }

    return {dname, args};
}

} // namespace

void DerivedQty::var_names(amrex::Vector<std::string>& plt_var_names)
{
    ioutils::add_var_names(plt_var_names, this->name(), this->num_comp());
}

DerivedQtyMgr::DerivedQtyMgr(const FieldRepo& repo) : m_repo(repo) {}

DerivedQty& DerivedQtyMgr::create(const std::string& key)
{
    auto qty_name = strip_spaces(key);

    // If this quantity is already registered return early
    if (contains(qty_name)) {
        return *m_derived_vec[m_obj_map[qty_name]];
    }

    auto tokens = parse_derived_qty(qty_name);
    m_derived_vec.emplace_back(
        DerivedQty::create(tokens.first, m_repo, tokens.second));
    m_obj_map[qty_name] = static_cast<int>(m_derived_vec.size()) - 1;

    return *m_derived_vec.back();
}

void DerivedQtyMgr::create(const amrex::Vector<std::string>& keys)
{
    for (const auto& qty : keys) {
        create(qty);
    }
}

void DerivedQtyMgr::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ALWAYS_ASSERT((scomp + num_comp()) <= fld.num_comp());

    int icomp = scomp;
    for (auto& qty : m_derived_vec) {
        (*qty)(fld, icomp);
        icomp += qty->num_comp();
    }
}

int DerivedQtyMgr::num_comp() const noexcept
{
    return std::accumulate(
        m_derived_vec.begin(), m_derived_vec.end(), 0,
        [](const int init, const std::unique_ptr<DerivedQty>& qty) {
            return init + qty->num_comp();
        });
}

bool DerivedQtyMgr::contains(const std::string& key) const noexcept
{
    auto it = m_obj_map.find(key);
    return (it != m_obj_map.end());
}

void DerivedQtyMgr::var_names(
    amrex::Vector<std::string>& plt_var_names) const noexcept
{
    for (const auto& qty : m_derived_vec) {
        qty->var_names(plt_var_names);
    }
}

} // namespace amr_wind
