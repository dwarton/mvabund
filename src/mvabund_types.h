// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// This file currently contains a copy of a recent/pending addition to Rcpp just so
// that we do NOT have to enforce a dependency on this very recent Rcpp because of it.
// In "due course", maybe in a year or two, the definition of 'Nullable' below can be
// removed when a corresponding dependency on Rcpp (> 0.12.0) is added.
//
// By keeing the same include guard, it is safe to have the file appear in both places
// in the meantime.
//
// We removed Doxygen tags which are in the original file as Rcpp is processed by Doxygen.
//
// The name of the file -- in the form PACKAGENAME_types.h -- is a convention used by
// Rcpp Attributes as we need the types defined here in the auto-generated interface.


// Nullable.h: Rcpp R/C++ interface class library -- SEXP container which can be NULL
//
// Copyright (C) 2015         Dirk Eddelbuettel
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#ifndef Rcpp_Nullable_h
#define Rcpp_Nullable_h

#include <RcppGSL.h>

namespace Rcpp {

    class Nullable {
    public:

        inline Nullable() : m_sexp(NULL), m_set(false) {}
        inline Nullable(SEXP t) : m_sexp(t), m_set(true) {}

        inline Nullable &operator=(SEXP sexp) {
            m_sexp = sexp;
            m_set = true;
            return *this;
        }

        inline operator SEXP() {
            checkIfSet();
            return m_sexp;
        }

        inline SEXP get() {
            checkIfSet();
            return m_sexp;
        }

        inline bool isNull() {
            checkIfSet();
            return Rf_isNull(m_sexp);
        }

        inline bool isNotNull() {
            return ! isNull();
        }

        inline bool isSet(void) { return m_set; }

    private:
        SEXP m_sexp;
        bool m_set;

        inline void checkIfSet(void) {
            if (!m_set) {
                throw ::Rcpp::exception("Not initialized");
            }
        }
    };
}

#endif
