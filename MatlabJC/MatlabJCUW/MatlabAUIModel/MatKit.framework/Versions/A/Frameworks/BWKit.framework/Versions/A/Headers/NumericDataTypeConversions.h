/*
  DataTypeConversions.h
  BWKit

  Created by Barry Wark on 2/22/08.
  Copyright 2008 Physion Consultants. All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 	- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 	- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 	- Neither the name of Physion Consultants nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#import <Cocoa/Cocoa.h>
#include <vector>
#include <memory>
#include <iterator>
#include <algorithm>
#include <functional>

using namespace std;

namespace bwkit {
    typedef char byte;
    
    /*!
     @functiongroup Useful functions
     */
    
    
    /*!
    @function
    @abstract   Swap the byte order of each element in a numeric array
    @discussion Swaps the endian byte order of each element. Obviously, the
     input should contain numeric data of the same type.
    @param      in NSData* containing numeric data of uniform type.
    @result     NSData* containing a copy of the input with all elements' endian
     order swapped.
    @templatefield T Type of the numeric data (e.g. double)
     */

    template<typename T>
    NSData *swap_numeric_data_byte_order(NSData *in); // usage: swap_numeric_data_byte_order<type>(data)
    
    /*!
     @function
     @abstract   Swap the byte order of each element in a numeric array in place.
     @discussion Swaps the endian byte order of each element in place. Obviously, the
     input should contain numeric data of the same type.
     @param      in NSMutableData* containing numeric data of uniform type. The data
     will be swapped in-place.
     @templatefield T Type of the numeric data (e.g. double)
     */
    template<typename T>
    void swap_numeric_data_byte_order(NSMutableData *in); //inplace swap for NSMutableData. usage: swap_numeric_data_byte_order<type>(data)
    
    
    /*!
     @function
     @abstract   Convert the type of elements in a numeric data array.
     @discussion Each element in the input is converted (via C type-casting) to
     the desired output type. No warning is produced at run time if information
     will be lost by the type conversion (though a compiler warning is likely).
     @param in NSData* containing numeric data of uniform type.
     @param inByteOrder Endian order of the input data.
     @param outByteOrder Desired endian order of the output data.
     @result     NSData* containing a copy of the input with all elements' type
     converted to OutType and endian order swapped if required to outByteOrder.
     @templatefield InType Type of the numeric data (e.g. double)
     @templatefield OutType Type of the desired output data (e.g. float)
     */
    template<typename InType, typename OutType>
    NSData *convert_numeric_data_type(NSData *in, CFByteOrder inByteOrder=NSHostByteOrder(), CFByteOrder outByteOrder=NSHostByteOrder()); //usage: newData = convert_numeric_data_type<InType,OutType>(inData, inByteOrder, outByteOrder)
    
    template<typename InputOutputIterator>
    void swap_byte_order(InputOutputIterator begin, 
                         InputOutputIterator end);
    
    
    /*!
     @functiongruop Less useful functions
     */
    
    
    template<typename T, typename U>
    auto_ptr<vector<U> > convert_data_type(auto_ptr<vector<T> > in);
    
    template<typename T>
    auto_ptr<vector<T> > swap_vector_byte_order(auto_ptr<vector<T> > vptr);
    
    template<typename T>
    T __byteswap(T v);
    
    /*!
     @functiongroup Converting std::vector <-> NSData*
     */
    template<typename T>
    auto_ptr<vector<T> > numeric_data_to_vector(NSData *d);
    
    template<typename T>
    NSData *vector_to_numeric_data(auto_ptr<vector<T> > v);
}


template<typename InType, typename OutType>
NSData *bwkit::convert_numeric_data_type(NSData *in, CFByteOrder inByteOrder=NSHostByteOrder(), CFByteOrder outByteOrder=NSHostByteOrder()) {
    auto_ptr<vector<OutType> > outPtr(convert_data_type<InType,OutType>(numeric_data_to_vector<InType>(in)));
    if(inByteOrder != outByteOrder) {
        swap_byte_order(outPtr->begin(), outPtr->end());
    }
    
    return vector_to_numeric_data(outPtr);
}

template<typename T>
NSData *bwkit::swap_numeric_data_byte_order(NSData *in) {
    using namespace bwkit;
    return vector_to_numeric_data(swap_vector_byte_order(numeric_data_to_vector<T>(in)));
}
    

template<typename T>
void bwkit::swap_numeric_data_byte_order(NSMutableData *in) {
    T *inPtr = (T*)[in mutableBytes];
    swap_byte_order(inPtr, inPtr+([in length]/sizeof(T)));
}

template<typename T>
auto_ptr<vector<T> > bwkit::swap_vector_byte_order(auto_ptr<vector<T> > vptr) {
    swap_byte_order(vptr->begin(), vptr->end());
    return vptr;
}

template<typename T, typename U>
auto_ptr<vector<U> > bwkit::convert_data_type(auto_ptr<vector<T> > in) {
    return auto_ptr<vector<U> >(new vector<U>(in->begin(), in->end()));
}

template<typename InputOutputIterator>
void bwkit::swap_byte_order(InputOutputIterator begin, 
                            InputOutputIterator end) {
    transform(begin, end, begin, pointer_to_unary_function<typename iterator_traits<InputOutputIterator>::value_type,
              typename iterator_traits<InputOutputIterator>::value_type> (bwkit::__byteswap));
}

template<typename T>
T bwkit::__byteswap(T v) {
    assert(sizeof(bwkit::byte)==1);
    bwkit::byte *vbytes = (bwkit::byte*)&v;
    reverse(vbytes, vbytes+sizeof(T));
    return v;
}

template<typename T>
auto_ptr<vector<T> > bwkit::numeric_data_to_vector(NSData *d) {
    auto_ptr<vector<T> > vptr(new vector<T>((T*)[d bytes], (T*)((T*)[d bytes]+([d length]/sizeof(T)))));
    
    return vptr;
}

template<typename T>
NSData *bwkit::vector_to_numeric_data(auto_ptr<vector<T> > vptr) {
    vector<T>& v = *vptr;
    NSData *result = [[NSData alloc] initWithBytes:&v[0] length:v.size()*sizeof(T)];
    
    return [result autorelease];
}
