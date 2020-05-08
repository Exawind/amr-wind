#ifndef IOMANAGER_H
#define IOMANAGER_H

#include <string>
#include <unordered_map>
#include <unordered_set>

#include "AMReX_Vector.H"

namespace amr_wind {

class CFDSim;
class Field;

class IOManager
{
public:
    IOManager(CFDSim&);

    ~IOManager();

    IOManager(const IOManager&) = delete;
    IOManager& operator=(const IOManager&) = delete;

    //! Process user inputs and perform one-time initialization tasks
    void initialize_io();

    //! Write all user-requested fields to disk
    void write_plot_file();

    //! Write all necessary fields for restart
    void write_checkpoint_file();

    //! Read all necessary fields for a restart
    void read_checkpoint_fields(const std::string& restart_file);

    //! Register a variable for output
    void register_output_var(const std::string& fname)
    {
        m_pltvars_default.insert(fname);
    }

    //! Register a variable for restart file
    void register_restart_var(const std::string& fname)
    {
        m_chkvars.insert(fname);
    }

    //! Convenience function to register a variable for both output and restart
    void register_io_var(const std::string& fname)
    {
        m_pltvars_default.insert(fname);
        m_chkvars.insert(fname);
    }

    const std::string& restart_file()
    { return m_restart_file; }

private:
    void write_header(const std::string&);

    void write_info_file(const std::string&);

    CFDSim& m_sim;

    //! Default output variables registered automatically in the code
    std::unordered_set<std::string> m_pltvars_default;

    //! Variables for output in checkpoint/restart files
    std::unordered_set<std::string> m_chkvars;

    //! Final list of fields to be output
    amrex::Vector<Field*> m_plt_fields;

    //! Final list of fields for restart
    amrex::Vector<Field*> m_chk_fields;

    //! Variable names (including components) for output
    amrex::Vector<std::string> m_plt_var_names;

    //! Prefix used for the plot file directories
    std::string m_plt_prefix{"plt"};

    //! Prefix used for the restart file directories
    std::string m_chk_prefix{"chk"};

    //! Restart file name
    std::string m_restart_file{""};

    //! Total number of variables (including components) output to plot file
    int m_plt_num_comp{0};

    //! Flag indicating whether default fields should be output
    bool m_output_default_vars{true};
};

} // namespace amr_wind

#endif /* IOMANAGER_H */
