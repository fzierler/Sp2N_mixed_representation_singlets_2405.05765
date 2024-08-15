function write_tex_table(name,data;insert_hline=[],no_header=false,extra_header=nothing)
    io = open(name,"w")
    rows, cols = size(data)
    # I have hardcoded centering of the table's contents
    table_layout = repeat("|c",cols)*"|"
    header = """\\begin{tabular}{$table_layout}
    \t\\hline
    """
    # either write a header and start with the second row or print the header
    write(io,header)
    isnothing(extra_header) || write(io,extra_header,"\n \\hline \\hline \n") 
    for i in 1:rows
        # if the gauge coupling changes insert two hlines for formatting
        if i ∈ insert_hline
            write(io,"\t\\hline\\hline\n")
        end
        # we could use some padding here for nicer formatting but for now
        # I only insert a '&' to make a minimal tex-compliant table
        write(io,"\t"*string(data[i,1])*"&")
        for j in 2:cols-1
            write(io,string(data[i,j])*"&")
        end
        write(io,string(data[i,cols])*"\\\\\n")
        #!no_header && i==1 && write(io,"\t\\hline\\hline\n")
    end
    write(io,"\t\\hline\\hline\n")
    write(io,"\\end{tabular}")
    close(io)
end
function write_tex_tables(tablepath,tex_tablepath)
    header_results   = L"   Label & $\beta$ & $N_t$ & $N_l$ & $am_0^{\rm f}$ & $am_0^{\rm as}$ & $am_{\eta^{\prime}_l}$ & $am_{\eta^{\prime}_h}$ & $am_{\rm PS}$ & $am_{\rm ps}$ & $am_{\rm V}$ & $am_{\rm v}$  \\\\"
    header_fitting   = L"   Label & $I_{\eta^{\prime}_l}$ & $I_{\eta^{\prime}_h}$ & $I_{\rm{PS}}$ & $I_{\rm{ps}}$ & $I_{\rm{V}}$ & $I_{\rm{v}}$ & $N_{\rm exp}$ & $\chi^2 / N_{\rm d.o.f.}$ & $\chi^2 / N_{\rm d.o.f.}$ & $\chi^2 / N_{\rm d.o.f.}$ & $\chi^2 / N_{\rm d.o.f.}$ & $\chi^2 / N_{\rm d.o.f.}$ & $\chi^2 / N_{\rm d.o.f.}$\\&&&&&&&&$ \eta^{\prime}_l$&$\eta^{\prime}_h$&${\rm PS}$&${\rm ps}$&$ {\rm V}$&$ {\rm v}$  \\\\"
    
    results = readdlm(joinpath(tablepath,"table_results.csv"),';',skipstart=1)
    fitting = readdlm(joinpath(tablepath,"table_fitting.csv"),';',skipstart=1)
    ispath(tex_tablepath) || mkpath(tex_tablepath)

    write_tex_table(joinpath(tex_tablepath,"table_results.tex"),results,extra_header=header_results)
    write_tex_table(joinpath(tex_tablepath,"table_fitting.tex"),fitting,extra_header=header_fitting)
end
function write_tex_table_t0_comparison(tablepath,tex_tablepath)
    header_results_l = L"   Label & $\beta$ & $N_t$ & $N_l$ & $am_0^{\rm f}$ & $am_0^{\rm as}$ & $am_{\eta^{\prime}_l} (t_0=1)$ & $am_{\eta^{\prime}_l} (t_0=2)$ & $am_{\eta^{\prime}_l} (t_0=3)$ & $am_{\eta^{\prime}_l} (t_0=4)$ & $am_{\eta^{\prime}_l} (t_0=5)$   \\\\"
    header_results_h = L"   Label & $\beta$ & $N_t$ & $N_l$ & $am_0^{\rm f}$ & $am_0^{\rm as}$ & $am_{\eta^{\prime}_h} (t_0=1)$ & $am_{\eta^{\prime}_h} (t_0=2)$ & $am_{\eta^{\prime}_h} (t_0=3)$ & $am_{\eta^{\prime}_h} (t_0=4)$ & $am_{\eta^{\prime}_h} (t_0=5)$   \\\\"
    
    results = readdlm(joinpath(tablepath,"table_results.csv"),';',skipstart=1)
    results_t0p1 = readdlm(joinpath(tablepath*"_t0p1","table_results.csv"),';',skipstart=1)
    results_t0p2 = readdlm(joinpath(tablepath*"_t0p2","table_results.csv"),';',skipstart=1)
    results_t0p3 = readdlm(joinpath(tablepath*"_t0p3","table_results.csv"),';',skipstart=1)
    results_t0p4 = readdlm(joinpath(tablepath*"_t0p4","table_results.csv"),';',skipstart=1)
    
    comparison_l = hcat(results[:,1:7],results_t0p1[:,7],results_t0p2[:,7],results_t0p3[:,7],results_t0p4[:,7])
    comparison_h = hcat(results[:,1:6],results[:,8],results_t0p1[:,8],results_t0p2[:,8],results_t0p3[:,8],results_t0p4[:,8])
    
    write_tex_table(joinpath(tex_tablepath,"table_results_t0_comparison_l.tex"),comparison_l,extra_header=header_results_l)
    write_tex_table(joinpath(tex_tablepath,"table_results_t0_comparison_h.tex"),comparison_h,extra_header=header_results_h)
end