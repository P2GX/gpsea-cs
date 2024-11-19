import io
import os
import typing
import matplotlib

def open_text_io_handle_for_reading(
    file: typing.Union[io.IOBase, str],
    encoding: str = "utf-8",
) -> io.TextIOWrapper:
    """
    Open a `io.TextIO` file handle based on `file`.

    :param file: a `str` or :class:`io.IOBase` to read from. If `str`, then the input is interpreted as a path
      to a local file.
      If `fh` is an IO wrapper, the function ensures we get a text wrapper that uses given encoding.
    :param encoding: encoding used to decode the input (`utf-8` by default).
    :return: the :class:`io.TextIOWrapper` wrapper.
    """
    if isinstance(file, str):
        # A path to local file
        return open(file, "r", encoding=encoding)
    elif isinstance(file, io.IOBase):
        if isinstance(file, (io.TextIOWrapper, io.TextIOBase)):
            return file
        elif isinstance(file, (io.BytesIO, io.BufferedIOBase)):
            return io.TextIOWrapper(file, encoding=encoding)

    raise ValueError(f"Unsupported type {type(file)}")


def open_text_io_handle_for_writing(
    file: typing.Union[io.IOBase, str],
    encoding: str = "utf-8",
) -> io.TextIOWrapper:
    """
    Open a `io.TextIO` file handle based on `file`.

    :param file: a `str` or :class:`io.IOBase` to write to. If `str`, then the input is interpreted as a path
      to a local file.
      If `fh` is an IO wrapper, the function ensures we get a text wrapper that uses given encoding.
    :param encoding: encoding used to encode the output (`utf-8` by default).
    :return: the :class:`io.TextIOWrapper` wrapper.
    """
    if isinstance(file, str):
        # A path to local file
        return open(file, "w", encoding=encoding)
    elif isinstance(file, io.IOBase):
        if isinstance(file, (io.TextIOWrapper, io.TextIOBase)):
            return file
        elif isinstance(file, (io.BytesIO, io.BufferedIOBase)):
            return io.TextIOWrapper(file, encoding=encoding)

    raise ValueError(f"Unsupported type {type(file)}")


def get_no_sig_fet_result_table(fet_d) -> typing.List[str]:
    general_info = fet_d.get("general_info") # dictionary of general results
    a_geno = general_info.get("a_genotype")
    b_geno = general_info.get("b_genotype")
    total_tests = general_info.get("n_tests_performed")
    latex_caption = f"""\
        Fisher Exact Test performed to compare HPO annotation frequency with respect to {a_geno} and {b_geno}."""
    interpr = fet_d.get("interpretation")
    if interpr is not None and len(interpr) > 1:
        latex_caption + interpr
    header = ["Genotype (A)", "Genotype (B)", "total tests performed", "significant results"]
    row = [a_geno, b_geno, str(total_tests), "0"]
    latex_rows = list()
    latex_rows.append("\\begin{subfigure}[b]{0.95\\textwidth}")
    latex_rows.append("\\centering")
    latex_rows.append("\\resizebox{\\textwidth}{!}{")
    latex_rows.append("\\begin{tabular}{llllrr}")
    latex_rows.append("\\toprule")
    latex_rows.append(" & ".join(header) + "\\\\")
    latex_rows.append("\\midrule")
    latex_rows.append(" & ".join(row) + "\\\\")    
    latex_rows.append("\\bottomrule")       
    latex_rows.append("\\end{tabular}")
    latex_rows.append("}")
    latex_rows.append("\\captionsetup{justification=raggedright,singlelinecheck=false}")
    latex_rows.append(f"\\caption{{ {latex_caption} }}")
    latex_rows.append("\\end{subfigure}")
    latex_rows.append("")
    latex_rows.append("\\vspace{2em}")
    latex_rows.append("")
    return latex_rows

def get_fet_result_table(fet_d) -> typing.List[str]:
    sig_result_list = fet_d.get("sig_result_list")
    if len(sig_result_list) == 0:
        return get_no_sig_fet_result_table(fet_d)
    general_info = fet_d.get("general_info") # dictionary of general results
    a_geno = general_info.get("a_genotype")
    b_geno = general_info.get("b_genotype")
    total_tests = general_info.get("n_tests_performed")
    latex_caption = f"""\
        Fisher Exact Test performed to compare HPO annotation frequency with respect to {a_geno} and {b_geno}. Total of
        {total_tests} tests were performed."""
    interpr = fet_d.get("interpretation")
    if interpr is not None and len(interpr) > 1:
        latex_caption + interpr
    header = ["HPO term", a_geno, b_geno, "p-value", "adj. p-value"]
    
    
    rows = list()
    for res in sig_result_list:
        with_geno_a = res.get("with_geno_a").replace("%", "\\%")
        with_geno_b = res.get("with_geno_b").replace("%", "\\%")
        row = [res.get("hpo_item"), with_geno_a, with_geno_b, format_for_latex(res.get("pval")), format_for_latex(res.get("adj_pval"))]
        rows.append(row)
    latex_rows = list()
    latex_rows.append("\\begin{subfigure}[b]{0.95\\textwidth}")
    latex_rows.append("\\centering")
    latex_rows.append("\\resizebox{\\textwidth}{!}{")
    latex_rows.append("\\begin{tabular}{llllrr}")
    latex_rows.append("\\toprule")
    latex_rows.append(" & ".join(header) + "\\\\")
    latex_rows.append("\\midrule")
    for r in rows:
        latex_rows.append(" & ".join(r) + "\\\\")    
    latex_rows.append("\\bottomrule")       
    latex_rows.append("\\end{tabular}")
    latex_rows.append("}")
    latex_rows.append("\\captionsetup{justification=raggedright,singlelinecheck=false}")
    latex_rows.append(f"\\caption{{ {latex_caption} }}")
    latex_rows.append("\\end{subfigure}")
    latex_rows.append("\\vspace{2em}")
    return latex_rows

def get_mono_result_table(mono_d) -> typing.List[str]:
    a_geno = mono_d.get("a_genotype")
    b_geno = mono_d.get("b_genotype")
    test_name = mono_d.get("test_name")
    var_name = mono_d.get("variable_name")
    desc = mono_d.get("description").replace("%", "\\%")
    pval = mono_d.get("pval")
    xrefs = mono_d.get("xrefs")
    interpr = mono_d.get("interpretation")
    latex_caption = f"{test_name} to compare {a_geno} and {b_geno} with respect to {var_name}."
    if interpr is not None and len(interpr) > 1:
        latex_caption = f"{latex_caption} {interpr}"
    header = ["Description", "Variable", "Genotype (A)", "Genotype (B)", "p-value", "xrefs"]
    row = [desc, var_name, a_geno, b_geno, format_for_latex(pval), xrefs]
    latex_rows = list()
    latex_rows.append("\\begin{subfigure}[b]{0.95\\textwidth}")
    latex_rows.append("\\captionsetup{justification=raggedright,singlelinecheck=false}")
    latex_rows.append("\\resizebox{\\textwidth}{!}{")
    latex_rows.append("\\begin{tabular}{llllrr}")
    latex_rows.append("\\toprule")
    latex_rows.append(" & ".join(header) + "\\\\")
    latex_rows.append("\\midrule")
    latex_rows.append(" & ".join(row) + "\\\\")    
    latex_rows.append("\\bottomrule")       
    latex_rows.append("\\end{tabular}")
    latex_rows.append("}")
   
    latex_rows.append(f"\\caption{{ {latex_caption} }}")
    latex_rows.append("\\end{subfigure}")
    latex_rows.append("")
    latex_rows.append("\\vspace{2em}")
    latex_rows.append("")
    return latex_rows

def format_for_latex(p_value):
    # Convert scientific notation to LaTeX format
    if "e" in p_value:
        coefficient, exponent = p_value.split('e')
        return f"${coefficient}\\times 10^{{{int(exponent)}}}$"
    else:
        return p_value


def process_latex_template(context_d: typing.Dict, 
                  mpt_fig: matplotlib.figure.Figure=None):

    output_dir = "../../supplement/img/"
    cohort = context_d.get("cohort_name").replace(" ", "_")
    cohort_name = context_d.get("cohort_name")
    interpretation = context_d.get("interpretation")
    xrefs = context_d.get("xrefs")
    lines = list()
    lines.append("\\begin{figure}[htbp]")
    lines.append("\\centering")
    if mpt_fig is not None:
        # save figure to output directory
        outname = f"{cohort}_protein_diagram-draft.pdf"
        output_file = os.path.join(output_dir, outname)
        if not os.path.isdir(output_dir):
            raise FileNotFoundError(f"Could not find output directory at {output_dir}")
        # Save the figure
        mpt_fig.tight_layout()
        mpt_fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        mpt_fig.savefig(output_file, format='pdf', bbox_inches='tight', pad_inches=0)
        print(f"Figure saved to {output_file}")
        latex_output_path = os.path.join("img", outname)
        # Subpanel 1: An image
        lines.append("\\begin{subfigure}[b]{0.95\\textwidth}")
        lines.append("\\centering")
        lines.append(f"\\includegraphics[width=\\textwidth]{{ {latex_output_path}}} ")
        lines.append("\\captionsetup{justification=raggedright,singlelinecheck=false}")
        lines.append(f"\\caption{{Distribution of variants in {cohort_name}}}")
        lines.append("\\end{subfigure}")
        lines.append("")
        lines.append("\\vspace{2em}")
        lines.append("")
    ## now show FET results, if any
    fet_results = context_d.get("fet_result_list")
    n_fet= context_d.get("n_fet_results")
    if n_fet > 0:
        for result in fet_results:
            lines.extend(get_fet_result_table(result))

    mono_results = context_d.get("mono_result_list")
    n_mono = context_d.get("n_mono_results")
    if n_mono > 0:
        for result in mono_results:
            lines.extend(get_mono_result_table(result))
    capt = context_d.get("caption")
    gene_capt = context_d.get("latex_gene_caption")
    lines.append(f"\\caption{{ {capt} {gene_capt}}}")
    lines.append("\\end{figure}")

    output_tex_file = f"../../supplement/tex/{cohort}_summary_draft.tex"
    fh = open(output_tex_file, "wt")
    for line in lines:
        fh.write(line + "\n")
    fh.close()
    print(f"Output to {output_tex_file}")
