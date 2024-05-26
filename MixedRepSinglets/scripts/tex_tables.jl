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