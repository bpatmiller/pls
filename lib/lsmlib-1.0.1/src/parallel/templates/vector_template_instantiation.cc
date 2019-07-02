/*
 * File:        vector_template_instantiation.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Explicit template instantiation for LSMLIB 
 */

#include <vector>

template class std::vector<bool>;
template class std::vector<char>;
template class std::vector<int>;
template class std::vector<float>;
template class std::vector<double>;
