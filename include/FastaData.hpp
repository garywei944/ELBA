// Created by Saliya Ekanayake on 1/7/19.

#ifndef LBL_DAL_FASTADATA_HPP
#define LBL_DAL_FASTADATA_HPP

#include <cstdlib>
#include <memory>
#include <vector>
#include <string>
#include "Types.hpp"
#include "TraceUtils.hpp"

/*!
 * Utility to store and retrive the data read from FASTA files.
 */
class FastaData {
public:
  /*!
   * Creates a new FastaData instance. Note, data could contain more
   * characters than required for this process. Therefore, l_start
   * and l_end identify the start and end (inclusive) of the data segment
   * relevant to this process. Note, l_end is expected to point to the last
   * character of the last sequence relevant to this process, so it doesn't
   * include a new line character at the end.
   *
   * While creating the FastData instance, any sequence with length less than
   * the k-mer size (k) is removed. Also, multi-line FASTA content is
   * modified in-place to be single lined. Moreover, any * characters that
   * can exist in the middle of a sequence is removed.
   *
   * The above removals will modify the l_end pointer as well the total number
   * of sequences used in the computation.
   * @param data Shared pointer to the raw character stream. Note, this may
   * contain more characters than necessary by the owning process.
   * @param l_start Starting offset of the data for this process.
   * @param l_end End offset (inclusive) of the data for this process.
   * @param k The k-mer size
   */
  FastaData(char *buff, ushort k, uint64_t l_start, uint64_t &l_end,
            const std::shared_ptr<TimePod> &tp, TraceUtils tu);

  /*!
   * Destructor for FastaData
   */
  ~FastaData();

  /*!
   * Returns a pointer to the raw FASTA data stream, and sets the
   * requested sequence's length, starting offset, and end offset (inclusive)
   * in the output parameters.
   * @param idx Requested sequence index local to this instance.
   * @param [out] len The length of the sequence.
   * @param [out] start_offset The start offset of the sequence.
   * @param [out] end_offset_inclusive The end offset (inclusive) of
   * the sequence.
   * @return A pointer to the raw character stream of the FASTA content.
   */
  char *get_sequence(uint64_t idx, ushort &len,
                     uint64_t &start_offset,
                     uint64_t &end_offset_inclusive);

  /*!
   * Returns a pointer to the sequence identifier in the raw FASTA data stream.
   * It also, sets the requested sequence identifier's length, starting offset,
   * end offset (inclusive) in the output parameters.
   *
   * @param idx Requested sequence index local to this instance.
   * @param [out] len The length of the sequence identifier.
   * @param [out] start_offset The start offset of the sequence identifier.
   * @param [out] end_offset_inclusive The end offset (inclusive) of the
   * sequence identifier
   * @return A pointer to the sequence identifier character stream.
   */
  char *get_sequence_id(uint64_t idx, ushort &len,
                        uint64_t &start_offset,
                        uint64_t &end_offset_inclusive);

  int get_sequence_len(uint64_t idx);

  std::string get_sequence_name(uint64_t idx);

  /*!
   * Sets the sub buffer size for sequences starting at
   * <tt>start_idx</tt> to and including <tt>end_idx_inclusive</tt>
   * in <tt>len</tt>. Also, sets the <tt>start_offset</tt>, and
   * <tt>end_offset_inclusive</tt> appropriately.
   * @param start_idx Index of the starting sequence for the interested range.
   * @param end_idx_inclusive Index (inclusive) of the ending sequence
   * in the range.
   * @param [out] len The lenght of the sub buffer.
   * @param [out] start_offset The start offset of the sequence range
   * (character offset) in the whole data array.
   * @param [out] end_offset_inclusive The end offset of the sequence range
   * (character offset) in the whole data array.
   */
  void buffer_size(uint64_t start_idx, uint64_t end_idx_inclusive,
                   uint64_t &len,
                   uint64_t &start_offset,
                   uint64_t &end_offset_inclusive);

  /*!
   *
   * @return The number of sequences local to this process.
   */
  uint64_t local_count();

  /*!
   * The original local sequence count including the number of sequences
   * removed due to length being less than k.
   * @return The total number of original sequences local to this instance.
   */
  uint64_t orig_local_count();

  /*!
   *
   * @return A vector containing the original local indices of the deleted
   * sequences.
   */
  uvec_64 *deleted_indices();

  /*!
   *
   * @return The (possibly modified) end offset of the data in the data buffer.
   * This could be different from the <tt>l_end</tt> parameter passed to the
   * constructor as some sequences might get pruned due to length being less
   * than <tt>k</tt>
   */
  uint64_t end_offset();

  /*!
   *
   * @return The start offset of the data in data buffer.
   */
  uint64_t start_offset();


  /*!
   *
   * @return a pointer to the data buffer.
   */
  const char *buffer();

  /*!
   * A helper method to print information of the FastaData instance.
   */
  void print();

private:
  std::shared_ptr<TimePod> tp = nullptr;
  TraceUtils tu;
  /*! Pointer to the raw character stream of FASTA content  */
  char *buff;
  /*! Number of sequences local to this process available in data stream */
  uint64_t l_seq_count;
  /*! Starting and ending offsets of the sequence data local to this instance
   * in the data stream
   */
  uint64_t l_start, l_end;
  /*! Offsets in the raw data stream to sequence ID starts */
  uvec_64 *id_starts = nullptr;
  /*! Offsets in the raw data stream to the sequence (actual bases) starts */
  uvec_64 *seq_starts = nullptr;

  /*! Deleted original local sequence indices.
   *  For example, if buff contains 10 sequences with original local
   *  indexes from [0,9] and indices 3, 4, and 7 are deleted then this would
   *  contain indices 3, 4, and 7. The local sequence indices will be from
   *  [0,6] and their mapping to original local sequence indices will be as
   *  (local idx --> original local index)
   *  0 --> 0, 1 --> 1, 2 --> 2,
   *  3 --> 5, 4 --> 6,
   *  5 --> 8, 6 --> 9
   *  To construct the original global index, one would need to add the original
   *  global sequence offset.
   */
  uvec_64 *del_idxs = nullptr;

  /*! The original local sequence count stored in this instance's buff.
   * This value is >= l_seq_count because we may remove sequences that are
   * less than k in length.
   */
  uint64_t orig_l_seq_count;
};


#endif //LBL_DAL_FASTADATA_HPP
