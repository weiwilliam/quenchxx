/*
 * (C) Copyright 2023 UCAR
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

namespace quenchxx {

// -------------------------------------------------------------------------------------------------

static std::map<std::string, std::vector<std::string>> vaderQuenchxxCustomCookbook() {
  return
    // pt: from t and pkz
    {{"potential_temperature",        {"AirPotentialTemperature_B"}},
    // P: from delp, from ps (and ak/bk)
    {"air_pressure_levels",          {"AirPressureAtInterface_B", "AirPressureAtInterface_A"}},
    // t: from p-pt and pt-base, from pt
    {"air_temperature",              {"AirTemperature_C", "AirTemperature_A"}},
    // p: from pe, from p-p and p-base
    {"air_pressure",                 {"AirPressure_A", "AirPressure_B"}},
    // rh:
    {"relative_humidity",            {"RelativeHumidity_A"}},
    // mr: from spfh
    {"humidity_mixing_ratio",        {"HumidityMixingRatio_A", "HumidityMixingRatio_B"}},
    // sulfmf:
    {"mass_fraction_of_sulfate_in_air",  {"SulfateMassFraction_A"}},
    // spfh: from mr
    {"specific_humidity",            {"HumidityMixingRatio_A"}},
    // ln(p) from pe
    {"ln_air_pressure_at_interface", {"LnAirPressureAtInterface_A"}},
    // qsat
    {"qsat",                         {"SaturationSpecificHumidity_A"}},
    // svp
    {"svp",                          {"SaturationVaporPressure_A"}},
    // dlsvpdT
    {"dlsvpdT",                      {"LogDerivativeSaturationVaporPressure_A"}},
    // p^kappa from pe and ln(p)
    {"air_pressure_to_kappa",        {"AirPressureToKappa_A"}},
    // delp: from p
    {"air_pressure_thickness",       {"AirPressureThickness_A"}},
    // pt: from t and ps
    {"potential_temperature",        {"AirPotentialTemperature_A"}},
    // ps: from delp
    {"surface_pressure",             {"SurfaceAirPressure_A"}},
    // tv: from t and q
    {"virtual_temperature",          {"AirVirtualTemperature_A"}}};
}

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
