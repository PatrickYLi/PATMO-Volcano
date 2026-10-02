import math

import patmo_string


def buildPhotoRates(network, options):
    all_rates = ""
    for reaction in network.photoReactions:
        all_rates += "!" + reaction.getVerbatim() + "\n"
        all_rates += "krate(:," + str(reaction.index) + ") = " + reaction.rate + "\n\n"

    zenith = float(options.zenith_angle)
    coef = float(options.TOA_para)
    if not math.isfinite(zenith) or not 0 <= zenith <= 180:
        raise ValueError("zenith_angle must be finite and between 0 and 180 degrees")
    if not math.isfinite(coef) or coef < 0:
        raise ValueError("TOA_para must be finite and nonnegative")
    mu = math.cos(math.radians(zenith)) if zenith < 90 else 0.0
    def double(value):
        return format(value, ".16e").replace("e", "d")
    patmo_string.fileReplaceBuild(
        "src_f90/patmo_photoRates.f90", "build/patmo_photoRates.f90",
        ["#PATMO_photoRates", "#PATMO_zenith_mu", "#PATMO_TOA_coef"],
        [all_rates, double(mu), double(coef)],
    )
