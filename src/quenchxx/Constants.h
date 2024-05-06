/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

// -------------------------------------------------------------------------------------------------

#include <string>
#include <vector>

// -------------------------------------------------------------------------------------------------

namespace quenchxx {

    // Interfaces for getting constants
    double getConstant(const std::string constName);
    std::vector<std::string> getAllConstantsNames();

}  // namespace quenchxx

// -------------------------------------------------------------------------------------------------
