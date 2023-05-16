/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt and Jan Christian Kässens,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of HybridGWAIS.
 *
 *    HybridGWAIS is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    HybridGWAIS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HybridGWAIS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef BUFFER_H
#define BUFFER_H

#include <vector>
#include <list>
#include <map>
#include <mutex>
#include <memory>

#ifdef USE_AD_FPGA
#include <admxrc3.h>
#include <adb3.h>
#endif

#include "BufferAllocator.h"

namespace hybridsys {

using namespace std;

template<typename T>
class BufferFactory;

template<typename T, template <typename> class Alloc>
class Buffer
{
    friend BufferFactory<Buffer<T,Alloc>>;

public:
    using value_type = T;
    using allocator = Alloc<T>;

#ifdef USE_AD_FPGA
    using BufferHandle = ADMXRC3_BUFFER_HANDLE;
    static constexpr BufferHandle BufferHandleInvalidValue = (unsigned)ADB3_HANDLE_INVALID_VALUE;
#else
    using BufferHandle = void*;
    static constexpr BufferHandle BufferHandleInvalidValue = nullptr;
#endif

    explicit Buffer(size_t buffer_size_, BufferHandle handle_ = BufferHandleInvalidValue)
        : buffer_size(buffer_size_), handle(handle_), content_length(0), metadata(NULL)
    {
        buffer.resize(buffer_size_);
    }

    ~Buffer() {}

    T *getData() {
        return buffer.data();
    }

    const T *getData() const {
        return buffer.data();
    }

    void   setContentLength(size_t length) { content_length = length; }
    size_t getContentLength() const { return content_length; }

    size_t getSize() const {
        return buffer_size;
    }

    BufferHandle getHandle() { return handle; }

    void setMetadata(shared_ptr<void> metadata_) { metadata = metadata_; }
    const shared_ptr<void> &getMetadata() const { return metadata; }

    // does not delete the contents, just sets the length to zero and releases the metadata
    void clear() {
    	setContentLength(0);
    	setMetadata(shared_ptr<void>());
    }

private:

    const size_t buffer_size;
    BufferHandle handle;
    vector<T, Alloc<T>> buffer;
    size_t content_length;
    shared_ptr<void> metadata;

};

using FPGABuffer = Buffer<char, PageAlignedAllocator>;
using CUDABuffer = Buffer<char, CUDAAllocator>;

}


#endif // BUFFER_H
