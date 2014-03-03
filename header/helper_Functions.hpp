//
//  helper_Functions.hpp
//  
//
//  Created by Sivaram Ambikasaran on 3/2/14.
//
//

#ifndef __helper_Functions_hpp__
#define __helper_Functions_hpp__

template <typename T>
T power(T x, unsigned n) {
        T y     =       1;
        for (int k=0; k<abs(n); ++k) {
                y*=x;
        }
        return y;
}

#endif /*(__helper_Functions_hpp__)*/