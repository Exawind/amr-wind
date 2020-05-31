#ifndef PDEFIELDS_H
#define PDEFIELDS_H

#include <iosfwd>

namespace amr_wind {

class FieldRepo;
class Field;

namespace pde {

struct PDEFields
{
    PDEFields(FieldRepo& repo_in, const std::string& var_name);

    FieldRepo& repo;
    Field& field;
    Field& mueff;

    Field& src_term;
    Field& diff_term;
    Field& conv_term;
};

}
}

#endif /* PDEFIELDS_H */
