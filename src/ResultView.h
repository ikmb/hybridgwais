/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt,
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

#ifndef RESULTVIEW_H
#define RESULTVIEW_H

#include <stdexcept>
#include <cinttypes>
#include <iterator>
#include <vector>
#include <ostream>

using namespace std;

class ResultDefaults {
public:
#ifdef SCORE_TYPE
    using default_score_type = SCORE_TYPE;
#else
    using default_score_type = double;
#endif
    using default_id_type = uint32_t;
};

// Remove __host__ and __device__ tags if we're not compiling CUDA code. This file will be used
// in CUDA and in non-CUDA compilation contexts. These definitions allow non-CUDA implementation files
// to include actual CUDA code used in here. Note that the actual CUDA code will not execute, this is
// just to make things compilable.
#ifndef __CUDACC__
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#ifndef __constant__
#define __constant__
#endif
#ifndef __global__
#define __global__
#endif
#endif

template<class I = ResultDefaults::default_id_type, class S = ResultDefaults::default_score_type> class Result;
template<class I = ResultDefaults::default_id_type, class S = ResultDefaults::default_score_type> class DeviceResult;
template<class I = ResultDefaults::default_id_type, class S = ResultDefaults::default_score_type> class ResultIterator;
template<class I = ResultDefaults::default_id_type, class S = ResultDefaults::default_score_type>
class ResultView
{
public:
    using id_type = I;
    using score_type = S;
    using iterator = ResultIterator<I,S>;
    using difference_type = ptrdiff_t;

    const size_t alignment = max(alignment_of<id_type>::value, alignment_of<score_type>::value);

    static constexpr size_t getAlignment() {
        return max(alignment_of<id_type>::value, alignment_of<score_type>::value);
    }

    ResultView(int idFields_, int scoreFields_, size_t bufferSize_, char *buffer_ = nullptr)
        : idFields(idFields_), scoreFields(scoreFields_), bufferSize(bufferSize_), buffer(buffer_) {
        auto id_block_size = idFields * sizeof(id_type);
        if(id_block_size % alignment != 0)
            id_block_size += (alignment - (id_block_size % alignment));

        auto score_block_size = scoreFields * sizeof(score_type);
        if(score_block_size % alignment != 0)
            score_block_size += (alignment - (score_block_size % alignment));

        resultSize = id_block_size + score_block_size;
    }

    ResultView(const ResultView&) = default; // default copy constructor

    __host__ __device__ size_t getResultSize() const { return resultSize; }
    size_t getResultCount() const { return bufferSize / resultSize; }

    __host__ __device__ size_t getScoreCount() const { return scoreFields; }
    __host__ __device__ size_t getIdCount() const { return idFields; }

    void setBuffer(char *new_buffer, size_t new_size) {
        this->buffer = new_buffer;
        this->bufferSize =new_size;
    }

    __host__ __device__ DeviceResult<I,S> getDeviceResult(size_t resultIndex) {
        ptrdiff_t id_block_size = idFields * sizeof(id_type);
        if(id_block_size % alignment != 0)
            id_block_size += (alignment - (id_block_size % alignment));

        char *target_buffer = buffer + resultIndex * resultSize;
        return DeviceResult<I,S>{ target_buffer, getIdCount(), id_block_size, getScoreCount() };
    }

    bool operator==(const ResultView<I,S> &rhs) const {
        return idFields == rhs.idFields
                && scoreFields == rhs.scoreFields
                && buffer == rhs.buffer
                && bufferSize == rhs.bufferSize;
    }

    ResultIterator<I,S> begin() { return ResultIterator<I,S>(0, *this); }
    ResultIterator<I,S> end() { return ResultIterator<I,S>(getResultCount(), *this); }

    friend class ResultIterator<I, S>;
private:
    unsigned idFields;
    unsigned scoreFields;
    size_t resultSize;
    size_t bufferSize;
    char *buffer;

};

template<class I, class S>
class DeviceResult {
public:
    using id_type = I;
    using score_type = S;

    __host__ __device__ DeviceResult(char *source, size_t idFields, ptrdiff_t scoreOffset, size_t scoreFields)
        : numIdFields(idFields), numScoreFields(scoreFields), scoreOfs(scoreOffset), bufferSource(source)
    {
    }

    __host__ __device__ void setID(size_t index, id_type new_id) {
       *reinterpret_cast<I*>(bufferSource + index * sizeof(id_type)) = new_id;
   }

   __host__ __device__ void setScore(size_t index, score_type new_score) {
       *reinterpret_cast<S*>(bufferSource + scoreOfs + index * sizeof(score_type)) = new_score;
   }

   __host__ __device__ unsigned getScoreFields() const { return numScoreFields; }
   __host__ __device__ unsigned getIdFields() const { return numIdFields; }

private:
   unsigned numIdFields;
   unsigned numScoreFields;
   ptrdiff_t scoreOfs;
   char *bufferSource;
};

template<class I, class S>
class Result {
public:
    using id_type = I;
    using score_type = S;

    Result(char *source, size_t idFields, ptrdiff_t scoreOffset, size_t scoreFields)
        : bufferSource(source), id(idFields), score(scoreFields)
    {
        for(unsigned i = 0; i < idFields; i++)
            id[i] = *reinterpret_cast<I*>(bufferSource + (i * sizeof(id_type)));
        source += scoreOffset;
        for(unsigned i = 0; i < scoreFields; i++)
            score[i] = *reinterpret_cast<S*>(source + (i * sizeof(score_type)));
    }

    const vector<id_type> &getID() const { return id;}
    const vector<score_type> &getScore() const { return score; }
    const id_type &getID(size_t index) const {
		if (index < id.size())
			return id[index];
		else
			throw out_of_range("Result definition does not have that much ID fields.");
	}
    const score_type &getScore(size_t index) const {
    	if (index < score.size())
    		return score[index];
    	else
			throw out_of_range("Result definition does not have that much score fields.");
    }
    void replaceID(size_t index, id_type id_) {
        if (index < id.size())
            id[index] = id_;
        else
            throw out_of_range("Result definition does not have that much ID fields");
    }

private:
    char *bufferSource;
    vector<id_type> id;
    vector<score_type> score;
};


/**
 * Compares two results component-wise. Result is undefined if both results have a different number of score values
 * or the comparison of their score types is undefined.
 */
template<class I, class S>
bool operator<(const Result<I, S> &lhs, const Result<I,S> &rhs) {
//    for(unsigned i = 0; i < lhs.getScore().size(); i++) {
        if(lhs.getScore(0) < rhs.getScore(0))
            return true;
        else if(lhs.getScore(0) <= rhs.getScore(0) && lhs.getScore(0) >= rhs.getScore(0)) { // equal scores
            vector<I> idsl = lhs.getID();
            vector<I> idsr = rhs.getID();
            for (size_t i = 0; i < idsl.size(); i++) {
                if (idsl[i] < idsr[i])
                    return true;
            }
        }
//    }
    return false; // both are equal if we reach this (read: not "less than")
}

template<class I, class S>
ostream& operator<<(ostream& out, const Result<I, S> &r) {
//    // write 1-based ID indices to output
//    for (const auto &id : r.getID()) {
//        out << (id+1) << "\t";
//    }

    // copy score fields to output
    // cannot use copy here because it also prints a trailing delimiter
    if (r.getScore().size())
        out << r.getScore(0);
    for(unsigned i = 1; i < r.getScore().size(); i++) {
    	out << "\t";
        out << r.getScore(i);
    }
    out << endl;
    return out;
}


template<class I, class S>
class ResultIterator
{
public:
    using id_type = I;
    using score_type = S;

    ResultIterator(ptrdiff_t offset_, ResultView<I,S> &view_)
        : offset(offset_), view(view_) {}

    const id_type &getID(size_t index) {
        if(index >= view.idFields)
            throw out_of_range("Result definition does not have that much ID fields");

        char *target_buf = view.buffer + view.getResultSize() * offset;

        return *reinterpret_cast<id_type*>(target_buf + (sizeof(id_type) * index));
    }

    const score_type &getScore(size_t index) {
        if(index >= view.scoreFields)
            throw out_of_range("Result definition does not have that much score fields");

        ptrdiff_t id_block_size = view.getIdCount() * sizeof(id_type);
        if(id_block_size % view.alignment != 0)
            id_block_size += (view.alignment - (id_block_size % view.alignment));

        char *target_buf = view.buffer + view.getResultSize() * offset;

        return *reinterpret_cast<score_type*>(target_buf + id_block_size + (sizeof(score_type) * index));
    }

    bool operator==(const ResultIterator<I,S> &rhs) const {
        return offset == rhs.offset && view == rhs.view;
    }

    bool operator!=(const ResultIterator<I,S> &rhs) const {
        return !(*this == rhs);
    }

    Result<I,S> operator*() const {
        ptrdiff_t id_block_size = view.getIdCount() * sizeof(id_type);
        if(id_block_size % view.alignment != 0)
            id_block_size += (view.alignment - (id_block_size % view.alignment));

        return Result<I,S>(view.buffer + (view.getResultSize() * offset), view.idFields, id_block_size, view.scoreFields);
    }

    ResultIterator& operator++() {
        ++offset;

        // also check underflow condition if offset is signed type
        if(offset < 0 || offset > static_cast<ptrdiff_t>(view.getResultCount()))
            throw out_of_range("Moving iterator past end of result list");
        return *this;
    }

    ResultIterator& operator--() {
        --offset;

        // also check overflow condition if offset wraps back to a positive range
        if(offset < 0 || offset > static_cast<ptrdiff_t>(view.getResultCount()))
            throw out_of_range("Moving iterator past end of result list");
        return *this;
    }

private:
    ptrdiff_t offset;
    ResultView<I,S> &view;
};

using ResultViewSpec = ResultView<>;

#endif // RESULTVIEW_H
