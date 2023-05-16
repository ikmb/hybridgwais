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

#ifndef BufferFactory_H
#define BufferFactory_H

#include <stack>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <functional>

#include "Buffer.h"
#include "BufferAllocator.h"

namespace hybridsys {

using namespace std;

class Hybridsys;
class FPGA;

template<typename T>
class BufferFactory
{

public:
    using buffer_type = T;

    explicit BufferFactory(size_t buffer_size_) : buffer_size(buffer_size_) {}

    ~BufferFactory() {
        while(!buffers.empty()) {
            delete buffers.top();
            buffers.pop();
        }
    }

    void preallocateBuffer() {
        buffers.push(new buffer_type(buffer_size));
//        buffer_count++;
    }

    BufferFactory(BufferFactory &&other) = default;
    BufferFactory& operator=(BufferFactory &&other) = default;

    BufferFactory(const BufferFactory &other) = delete;
    BufferFactory& operator=(const BufferFactory &other) = delete;

    shared_ptr<buffer_type> getIfMoreThanXAvailable(size_t x, bool wait = true) {

        shared_ptr<buffer_type> ret;

        unique_lock<mutex> lock(this->stack_guard);

        if(buffers.size() <= x && !wait) {
            // sorry, not enough buffers available but you insisted not to wait
            lock.unlock();

        } else if( (buffers.size() <= x && wait) || buffers.size() > x ) {

        	// wait for a buffer to become available, if necessary
            if(buffers.size() <= x) {
				auto t_start = chrono::high_resolution_clock::now();
				stack_cv.wait(lock, [&](){ return buffers.size() > x; });
				auto t_end = chrono::high_resolution_clock::now();
				chrono::milliseconds dur = chrono::duration_cast<chrono::milliseconds>(t_end - t_start);
				buffer_wait_time += dur;
            }

            // Deleter takes only one argument, this->put() takes two (this and the buffer pointer),
            // so bind the first one (this) to get a 1-ary function.
            ret = shared_ptr<buffer_type>(buffers.top(),
                                                   bind(&BufferFactory::put, this, placeholders::_1));
            buffers.pop();

            lock.unlock();
        }

        return ret;

    }

    shared_ptr<buffer_type> get(bool wait = true) {
        return getIfMoreThanXAvailable(0, wait);
    }


    size_t getBufferSize() const { return buffer_size; }
    //size_t getBufferCount() const { return buffer_count; }

private:

    size_t buffer_size;
//    size_t buffer_count = 0;
    stack<buffer_type*> buffers;
    mutex stack_guard;
    condition_variable stack_cv;
    chrono::milliseconds buffer_wait_time;

    void put(buffer_type *buffer) {
        lock_guard<mutex> guard(this->stack_guard);

        if(buffer != NULL && buffer != nullptr) {
            buffer->clear();
            buffers.push(buffer);
     //       cerr << "[put " << buffers.size() << "]" << flush;
            stack_cv.notify_one();
        } else {
     //       cerr << "[put null]" << flush;
        }
    }



};

}

#endif // BufferFactory_H
