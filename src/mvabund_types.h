// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// This file contains a copy of this recent addition to Rcpp just so that we do NOT
// have to enforce a dependency on this very recent Rcpp because of it. In "due
// course", maybe in a year or two, the definition of 'Nullable' below can be
// removed when a corresponding dependency on Rcpp (>= 0.12.1) is added.
//
// By keeing the same include guard, it is safe to have the file appear in both places
// in the meantime.
//
// The name of the file -- in the form PACKAGENAME_types.h -- is a convention used by
// Rcpp Attributes as we need the types defined here in the auto-generated interface.

#ifndef Rcpp_Nullable_h
#define Rcpp_Nullable_h

// We also need to (re-)include Rcpp.h to define traits etc. The remainder 
// of the file is unchanged relative to Rcpp's inst/include/Rcpp/Nullable.h
#include <Rcpp.h>

// We looked into the safe_bool_idiom [1] but found that more trouble than is
// warranted here.  We first and foremost want an operator SEXP() which got in
// the way of redefining operator bool.
// [1] http://www.artima.com/cppsource/safebool.html)


namespace Rcpp {

    template<class T>
    class Nullable {
    private:
        template<class U>
        friend class InputParameter;

        template<class U>
        friend class traits::Exporter;

    public:

        /**
         * Empty no-argument constructor of a Nullable object
         *
         * Assigns (R's) NULL value, and sets validator to FALSE
         */
        inline Nullable() : m_sexp(R_NilValue), m_set(false) {}

        /**
         * Template constructor of a Nullable object
         *
         * Assigns object, and set validator to TRUE
         */

        inline Nullable(const T &t) : m_sexp(t),  m_set(true) {}

    protected:

        /**
         * Standard constructor of a Nullable object
         *
         * @param SEXP is stored
         */
        inline Nullable(SEXP t) {
            m_sexp = t;
            m_set = true;
        }

    public:

        /**
         * Copy constructor for Nullable object
         *
         * @param SEXP is used to update internal copy
         */
        inline Nullable &operator=(SEXP sexp) {
            m_sexp = sexp;
            m_set = true;
            return *this;
        }

        /**
         * operator SEXP() to return nullable object
         *
         * @throw 'not initialized' if object has not been set
         */
        inline operator SEXP() {
            checkIfSet();
            return m_sexp;
        }

        /**
         * get() accessor for object
         *
         * @throw 'not initialized' if object has not been set
         */
        inline SEXP get() {
            checkIfSet();
            return m_sexp;
        }

        /**
         * Boolean test for NULL
         *
         * @throw 'not initialized' if object has not been set
         */
        inline bool isNull() const {
            checkIfSet();
            return Rf_isNull(m_sexp);
        }

        /**
         * Boolean test for not NULL
         *
         * @throw 'not initialized' if object has not been set
         */
        inline bool isNotNull() const {
            return ! isNull();
        }

        /**
         * Test function to check if object has been initialized
         *
         */
        inline bool isSet(void) const { return m_set; }

    private:
        SEXP m_sexp;
        bool m_set;

        inline void checkIfSet(void) const {
            if (!m_set) {
                throw ::Rcpp::exception("Not initialized");
            }
        }
    };
}

#endif
