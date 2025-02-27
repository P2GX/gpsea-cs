


def format_p(number) -> str: 
    n = float(number)
    if n > 0.001:
        return number
    else:
        return f"${n:.1e}".replace('e', r' \times 10^{') + "}$"