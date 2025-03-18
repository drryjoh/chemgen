import math

def is_integer_to_precision(x, tol=1e-9):
    """
    Check if the floating-point number x is an integer to machine precision.
    
    Parameters:
    - x: The floating-point number to check.
    - tol: The tolerance for checking. Default is 1e-9, which is suitable for many cases.
    
    Returns:
    - True if x is close to an integer, False otherwise.
    """
    return math.isclose(x, round(x), abs_tol=tol)

def raise_to_power(raise_to, power):
    if is_integer_to_precision(power, 1e-6):
        power_integer = int(power)
        if power_integer == 1:
            return raise_to
        else:
            return f"pow_gen{power_integer}({raise_to})"
    else:
        return f"pow_gen({raise_to}, {power})"

def draise_to_power_chain(raise_to, draise_to, power):
    if is_integer_to_precision(power, 1e-6):
        power_integer = int(power)
        if power_integer == 1:
            return draise_to
        else:
            return f"multiply(dpow_gen{power_integer}_da({raise_to}), {draise_to})"
    else:
        return f"multiply(dpow_gen_da({raise_to}, {power}), {draise_to})"

def draise_to_power(raise_to, power):
    if is_integer_to_precision(power, 1e-6):
        power_integer = int(power)
        if power_integer == 1:
            return "1"
        else:
            return f"dpow_gen{power_integer}_da({raise_to})"
    else:
        return f"dpow_gen_da({raise_to}, {power})"