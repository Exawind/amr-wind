#ifndef CONSOLE_IO_H
#define CONSOLE_IO_H

#include <iostream>
#include "AMReX_MLMG.H"

namespace amr_wind {
namespace io {

void print_banner(MPI_Comm /*comm*/, std::ostream& /*out*/);

void print_usage(MPI_Comm /*comm*/, std::ostream& /*out*/);

void print_error(MPI_Comm /*comm*/, const std::string& msg);

void print_summary(std::ostream&);

void print_mlmg_header(const std::string& /*key*/);

void print_mlmg_info(const std::string& solve_name, const amrex::MLMG& mlmg);

void print_tpls(std::ostream& /*out*/);

} // namespace io
} // namespace amr_wind

#endif /* CONSOLE_IO_H */
