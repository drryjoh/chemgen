def arrhenius_text(i, A, B, E, configuration = None):
    if configuration == None:
        return f"auto reaction_{i}(const auto& temperature) const {{ return arrhenius({A}, {B}, {E}, temperature)}}"
    else:
        return_text = "{device_option}\n{scalar_function} call_forward_reaction_{i}({scalar_parameter} temperature) {const_option} {{ return arrhenius({scalar_cast}({A}), {scalar_cast}({B}), {scalar_cast}({E}), temperature);}}"
        return return_text.format(**vars(configuration), i = i, A = A, E = E, B = B)

def darrhenius_text(i, A, B, E, configuration = None):
    if configuration == None:
        return f"auto dreaction_dtemperature_{i}(const double& temperature) const {{ return darrhenius_dtemperature({A}, {B}, {E}, temperature)}}"
