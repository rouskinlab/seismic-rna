// Required for Python <=3.12, before #include <Python.h>
#define PY_SSIZE_T_CLEAN
// Required for integration with Python
#include <Python.h>

// Standard C library
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// Exceptions

static PyObject *RelateError;


// Numbers

static const int DECIMAL = 10;


// Base encodings

static const char BASEA = 'A';
static const char BASEC = 'C';
static const char BASEG = 'G';
static const char BASET = 'T';


// Relationship encodings

static const unsigned char MATCH = '\x01';
static const unsigned char DELET = '\x02';
static const unsigned char INS_5 = '\x04';
static const unsigned char INS_3 = '\x08';
static const unsigned char SUB_A = '\x10';
static const unsigned char SUB_C = '\x20';
static const unsigned char SUB_G = '\x40';
static const unsigned char SUB_T = '\x80';
static const unsigned char SUB_N = SUB_A | SUB_C | SUB_G | SUB_T;
static const unsigned char ANY_N = SUB_N | MATCH;


// SAM file parsing

static const char *SAM_SEP = "\t\n";
static const char CIG_ALIGN = 'M';
static const char CIG_DELET = 'D';
static const char CIG_INSRT = 'I';
static const char CIG_MATCH = '=';
static const char CIG_SCLIP = 'S';
static const char CIG_SUBST = 'X';
static const unsigned long FLAG_PAIRED = 1;
static const unsigned long FLAG_REV = 16;
static const unsigned long FLAG_1ST = 64;
static const unsigned long FLAG_2ND = 128;
static const unsigned long MAX_FLAG = 4095;  // 4095 = 2^12 - 1


/* Return the smallest of two arguments. */
static inline size_t min(size_t a, size_t b)
{
    return (a < b) ? a : b;
}


/* Return the largest of two arguments. */
static inline size_t max(size_t a, size_t b)
{
    return (a > b) ? a : b;
}


/*
Parse a string (str) and store the value in an integer (number).
A basic wrapper around string-to-unsigned-long (strtoul) that also
determines whether a return value of 0 is correct (e.g. if str == "0")
or if it is because parsing failed (strtoul returns 0 upon failure).
Overflows are NOT checked for because in practice none of the values in
the SAM file will be so large as to overflow an unsigned long type.

Parameters
----------
ulong
    Non-NULL pointer to number in which to store the result
str
    Non-NULL pointer to string from which to parse the number

Returns
-------
int
    0 if successful, otherwise an error code (> 0)
*/
static int parse_ulong(unsigned long *ulong, const char *str)
{
    // Neither ulong nor str may be NULL.
    assert(ulong != NULL);
    assert(str != NULL);
    // Parse str as a base-10 integer; store its numeric value in ulong.
    char *endptr;
    *ulong = strtoul(str, &endptr, DECIMAL);
    // Check if the parse failed to process any characters.
    if (endptr == str)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse unsigned integer");
        return 1;
    }
    return 0;
}


typedef struct
{
    // Parsed attributes
    const char *name;    // read name
    unsigned long flag;  // bitwise flag
    const char *ref;     // reference name
    size_t pos;          // mapping position in the ref (1-indexed)
    unsigned long mapq;  // mapping quality
    const char *cigar;   // CIGAR string
    const char *seq;     // read sequence
    const char *qual;    // read quality
    // Calculated attributes
    const char *end;     // pointer to address just after 3' end of read
    size_t len;          // read length
    int paired;          // whether read is paired
    int reverse;         // whether read is reverse complemented
    int first;           // whether read is the 1st read
    int second;          // whether read is the 2nd read
    size_t ref_end5;     // 5' end of the read in the ref (1-indexed)
    size_t ref_end3;     // 3' end of the read in the ref (1-indexed)
} SamRead;


/*
Parse one line from a SAM file and store the information in each field
in one of the field arguments passed to this function.

Parameters
----------
read
    Non-nullable pointer to the SAM read struct in which to store the
    parsed information.
line
    Non-nullable pointer to the text of a SAM-formatted line to parse.
    It must be tab-delimited and contain at least 11 fields (below).

Returns
-------
int
    0 if successful, otherwise an error code (> 0)
*/
static int parse_sam_line(SamRead *read, const char *line)
{
    assert(read != NULL);
    assert(line != NULL);
    char *end = (char *)line;  // End of the current field.
    unsigned long temp_ulong;  // Hold parsed numbers before casting.

    // Read name
    if ((read->name = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse read name");
        return 1;
    }

    // Bitwise flag
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end)))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse SAM flag");
        return 1;
    }
    read->flag = temp_ulong;
    if (read->flag > MAX_FLAG)
    {
        PyErr_SetString(PyExc_ValueError, "SAM flag is too large");
        return 1;
    }
    // Individual flag bits
    read->paired = (read->flag & FLAG_PAIRED) > 0;
    read->reverse = (read->flag & FLAG_REV) > 0;
    read->first = (read->flag & FLAG_1ST) > 0;
    read->second = (read->flag & FLAG_2ND) > 0;

    // Reference name
    if ((read->ref = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse reference name");
        return 1;
    }

    // Mapping position
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end)))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse mapping position");
        return 1;
    }
    read->pos = (size_t)temp_ulong;

    // Mapping quality
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end)))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse mapping quality");
        return 1;
    }
    read->mapq = temp_ulong;

    // CIGAR string
    if ((read->cigar = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse CIGAR string");
        return 1;
    }

    // Next reference (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse next reference name");
        return 1;
    }

    // Next position (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse next mapping position");
        return 1;
    }

    // Template length (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse template length");
        return 1;
    }

    // Read sequence
    if ((read->seq = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse read sequence");
        return 1;
    }

    // Read end position
    // Subtract 1 because end points to one after the \t delimiter, but
    // read->end must point to the tab delimiter itself, i.e.
    // ACGT FFFF      (length = 4)
    // ^read->seq
    //     ^read->end (= 4 + read->seq)
    //      ^end      (= 5 + read->seq)
    read->end = end - 1;

    // Read length (using pointer arithmetic)
    read->len = read->end - read->seq;

    // Read quality
    if ((read->qual = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse read quality");
        return 1;
    }

    // Lengths of read and quality strings must match.
    // Subtract 1 from end for the same reason as in "Read end position"
    if (end != NULL && read->len != (size_t)((end - 1) - read->qual))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Sequence and quality strings differ in length");
        return 1;
    }

    // Initialize the 5' and 3' ends to placeholder 0 values.
    read->ref_end5 = 0;
    read->ref_end3 = 0;

    // Parsing the line succeeded.
    return 0;
}


static int validate_read(const SamRead *read,
                         const char *ref,
                         unsigned long min_mapq)
{
    assert(read != NULL);
    assert(read->ref != NULL);
    assert(ref != NULL);
    if (read->mapq < min_mapq)
    {
        PyErr_SetString(PyExc_ValueError, "Mapping quality is insufficient");
        return 1;
    }
    if (strcmp(read->ref, ref))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Reference name does not match name of SAM file");
        return 1;
    }
    return 0;
}


static int validate_pair(const SamRead *read1,
                         const SamRead *read2)
{
    assert(read1 != NULL);
    assert(read2 != NULL);
    assert(read1->name != NULL);
    assert(read2->name != NULL);
    if (strcmp(read1->name, read2->name))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 have different reference names");
        return 1;
    }
    if (!(read1->paired && read2->paired))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 are not both paired-end");
        return 1;
    }
    if (!(read1->first && read2->second))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 are not marked as 1st and 2nd");
        return 1;
    }
    if (read1->reverse == read2->reverse)
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 aligned in the same orientation");
        return 1;
    }
    return 0;
}


/* Encode a base character as a substitution. */
static unsigned char encode_subs(char base)
{
    switch (base)
    {
        case BASET:
            return SUB_T;
        case BASEG:
            return SUB_G;
        case BASEC:
            return SUB_C;
        case BASEA:
            return SUB_A;
        default:
            return SUB_N;
    }
}


/*
Encode a match as a match byte (if the read quality is sufficient),
otherwise as an ambiguous byte (otherwise).
*/
static unsigned char encode_match(char read_base,
                                  char read_qual,
                                  unsigned char min_qual)
{
    return (read_qual >= min_qual) ? MATCH : (ANY_N ^ encode_subs(read_base));
}


static unsigned char encode_relate(char ref_base,
                                   char read_base,
                                   char read_qual,
                                   unsigned char min_qual)
{
    if (read_qual < min_qual) {return ANY_N ^ encode_subs(ref_base);}
    return (ref_base == read_base) ? MATCH : encode_subs(read_base);
}


// CIGAR string operations

typedef struct
{
    const char *op;    // type of operation
    size_t len;  // length of operation
} CigarOp;


/*
Find the next operation in a CIGAR string. Store the kind of operation
in cigar->op and the length of the operation in cigar->len. If parsing
does not return a valid operation (which must have cigar->len > 0 and
cigar->op != NULL), then set cigar->op to NULL. This event happens when
the end of the CIGAR string is reached and does not signify an error.

Parameters
----------
op
    Non-nullable pointer to a CIGAR operation struct in which to store
    the parsed information from the CIGAR string

Returns
-------
int
    0 if successful, otherwise an error code >0
*/
static int get_next_cigar_op(CigarOp *cigar)
{
    // Parse as many characters of text as possible to a number, cast
    // to uint32, and store in cigar->len (length of CIGAR operation).
    // Parsing must start one character ahead of the last character to
    // be read from the CIGAR string (which is cigar->op, the type of
    // the previous operation); thus, parsing starts at cigar->op + 1.
    // The character after the last digit to be parsed is the type of
    // the new operation; thus, after cigar->len is parsed, cigar->op
    // is finally pointed at the character after the last parsed digit.
    cigar->len = (size_t)strtoul(cigar->op + 1, (char **)&cigar->op, DECIMAL);
    if (cigar->len == 0)
    {
        // The parse returns a length of 0 if it fails. This is normal
        // when the end of the string is reached, but should not happen
        // before then. If the character to which the CIGAR operation
        // points evaluates to true, then it is not the null terminator
        // of a string, so the end of end of the CIGAR string has not
        // yet been reached, which is an error.
        if (*cigar->op) {return 402;}
        // Otherwise, simply point the operation to NULL to signal the
        // normal end of the CIGAR string.
        cigar->op = NULL;
    }
    return 0;
}


// Insertions and deletions

static const int DELETION = 0;
static const int INSERTION = 1;
static const size_t MAX_NUM_PODS = 16;
static const size_t MAX_POD_SIZE = 64;


static unsigned char get_ins_rel(int insert3)
{
    return insert3 ? INS_3 : INS_5;
}


static size_t calc_lateral5(size_t lateral3)
{
    return (lateral3 > 0) ? lateral3 - 1 : 0;
}


static size_t get_lateral(size_t lateral3, int insert3)
{
    return insert3 ? lateral3 : calc_lateral5(lateral3);
}


typedef struct
{
    int insert;
    size_t opposite;
    size_t lateral3;
} Indel;


typedef struct
{
    size_t n;
    int insert;
    Indel indels[MAX_POD_SIZE];
    size_t size;
} IndelPod;


typedef struct
{
    IndelPod pods[MAX_NUM_PODS];
    size_t size;
} IndelPodArray;


static int add_indel(IndelPodArray *pods,
                     int insert,
                     size_t opposite,
                     size_t lateral3)
{
    IndelPod *pod = NULL;
    if ((pods->size == 0) || (pods->pods[pods->size - 1].insert != insert))
    {
        // Check if the maximum number of pods will be exceeded.
        if (pods->size >= MAX_NUM_PODS)
        {
            PyErr_SetString(PyExc_ValueError,
                            "Read has too many pods of indels");
            return 1;
        }
        pod = &(pods->pods[pods->size]);
        // Initialize the members of a new pod.
        pod->n = pods->size;
        pod->insert = insert;
        pod->size = 0;
        // Increment the number of pods.
        pods->size++;
    }
    // Choose the last pod.
    pod = &(pods->pods[pods->size - 1]);
    // Check if the maximum number of indels in a pod will be exceeded.
    if (pod->size >= MAX_POD_SIZE)
    {
        PyErr_SetString(PyExc_ValueError,
                        "Too many indels in one pod");
        return 1;
    }
    Indel *indel = &(pod->indels[pod->size]);
    // Initialize the members of a new indel.
    indel->insert = insert;
    indel->opposite = opposite;
    indel->lateral3 = lateral3;
    // Increment the number of indels in the pod.
    pod->size++;
    return 0;
}


/*
Compute the relationships of a SamRead.

Parameters
----------
read
    Non-nullable pointer to the SamRead.
ref_seq
    Non-nullable pointer to the sequence of the reference within the
    region of interest. The sequence may contain only the characters
    'A', 'C', 'G', and 'T' (lowercase not allowed), and its length must
    equal sect_len, otherwise the behavior is undefined.
sect_len
    Length of the section. Must equal lengths of muts and sect_seq,
    otherwise the behavior is undefined and may cause memory violations.
sect_end5
    Position of the 5' end of the section with respect to the beginning
    of the entire reference sequence (1-indexed). Must be positive.
read
    Read from a SAM file.
min_qual
    Minimum ASCII-encoded quality score to accept a base call.
ambid
    Whether to compute and label ambiguous insertions and deletions.

Returns
-------
On success: 0
On failure: >0
*/
static int calc_rels_read(unsigned char *rels,
                          SamRead *read,
                          const char *ref_seq,
                          size_t ref_len,
                          unsigned char min_qual,
                          int insert3,
                          int ambindel,
                          size_t clip_end5,
                          size_t clip_end3)
{
    // Validate the arguments.
    int error;
    assert(rels != NULL);
    assert(read != NULL);
    assert(ref_seq != NULL);
    if (read->len == 0)
    {
        PyErr_SetString(PyExc_ValueError, "Length of read sequence is 0");
        return 1;
    }
    if (ref_len == 0)
    {
        PyErr_SetString(PyExc_ValueError, "Length of reference sequence is 0");
        return 1;
    }
    if (read->pos == 0)
    {
        PyErr_SetString(PyExc_ValueError, "Mapping position is 0");
        return 1;
    }
    if (read->pos > ref_len)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "Mapping position is greater than length of reference sequence"
        );
        return 1;
    }

    // Positions in the reference and read (0-indexed).
    size_t ref_pos = read->pos - 1;
    size_t read_pos = 0;
    // 5' and 3' ends of the read, in the read coordinates (1-indexed).
    size_t read_end5 = 1;
    size_t read_end3 = read->len;
    // CIGAR operation stopping point.
    size_t cigar_op_stop_pos;

    // Positions of deletions and insertions.
    IndelPodArray pods;
    pods.size = 0;
    
    // Initialize the CIGAR operation so that it points just before the
    // CIGAR string.
    CigarOp cigar;
    cigar.op = read->cigar - 1;
    // Read the first operation from the CIGAR string; catch errors.
    if (get_next_cigar_op(&cigar)) {return 402;}
    // Return an error if there were no operations in the CIGAR string.
    if (cigar.op == NULL) {return 401;}

    // Read the entire CIGAR string one operation at a time.
    while (cigar.op != NULL)
    {
        // Decide what to do based on the current CIGAR operation.
        switch (*cigar.op)
        {
    
        case CIG_MATCH:
            // The read and reference match over the entire operation.
            cigar_op_stop_pos = ref_pos + cigar.len;
            if (cigar_op_stop_pos > ref_len) {return 404;}
            if (read_pos + cigar.len > read->len) {return 404;}
            while (ref_pos < cigar_op_stop_pos)
            {
                rels[ref_pos] = encode_match(read->seq[read_pos],
                                             read->qual[read_pos],
                                             min_qual);
                ref_pos++;
                read_pos++;
            }
            break;
        
        case CIG_ALIGN:
        case CIG_SUBST:
            // The read and reference have matches or substitutions over
            // the entire operation.
            cigar_op_stop_pos = ref_pos + cigar.len;
            if (cigar_op_stop_pos > ref_len) {return 404;}
            if (read_pos + cigar.len > read->len) {return 404;}
            while (ref_pos < cigar_op_stop_pos)
            {
                rels[ref_pos] = encode_relate(ref_seq[ref_pos],
                                              read->seq[read_pos],
                                              read->qual[read_pos],
                                              min_qual);
                ref_pos++;
                read_pos++;
            }
            break;
        
        case CIG_DELET:
            // The portion of the reference sequence corresponding to
            // the operation is deleted from the read.
            cigar_op_stop_pos = ref_pos + cigar.len;
            if (cigar_op_stop_pos > ref_len) {return 404;}
            if (read_pos == 0 || read_pos >= read->len - 1) {return 1;}
            while (ref_pos < cigar_op_stop_pos)
            {
                rels[ref_pos] = DELET;
                error = add_indel(&pods, DELETION, ref_pos, read_pos);
                if (error) {return error;}
                ref_pos++;
            }
            break;
        
        case CIG_INSRT:
            // The read contains an insertion of one or more bases that
            // are not present in the reference sequence. 
            cigar_op_stop_pos = read_pos + cigar.len;
            if (cigar_op_stop_pos > read->len) {return 404;}
            if (ref_pos == 0 || ref_pos >= ref_len - 1) {return 1;}
            while (read_pos < cigar_op_stop_pos)
            {
                error = add_indel(&pods, INSERTION, read_pos, ref_pos);
                if (error) {return error;}
                read_pos++;
            }
            break;
        
        case CIG_SCLIP:
            // Bases were soft-clipped from the 5' or 3' end of the read
            // during alignment. Like insertions, they consume the read
            // but not the reference.
            cigar_op_stop_pos = read_pos + cigar.len;
            if (cigar_op_stop_pos > read->len) {return 404;}
            if (read_pos == 0)
            {
                // This is the soft clip from the 5' end of the read.
                if (read_end5 != 1) {return 1;}
                read_end5 += cigar.len;
            }
            else
            {
                // This is the soft clip from the 3' end of the read.
                if (read_end3 != read->len) {return 1;}
                read_end3 -= cigar.len;
            }
            read_pos = cigar_op_stop_pos;
            break;
        
        default:
            // The CIGAR operation was not recognized.
            return 403;
        }

        // Read the next operation from the CIGAR string; catch errors.
        if (get_next_cigar_op(&cigar)) {return 402;}
    }

    // Verify that the sum of all CIGAR operations that consumed the
    // read equals the length of the read. The former equals read_pos
    // because for each CIGAR operation that consumed the read, the
    // length of the operation was added to read_pos.
    if (read_pos != read->len) {return 405;}
    // Add insertions to rels.
    const unsigned char ins_rel = get_ins_rel(insert3);
    for (size_t p = 0; p < pods.size; p++)
    {
        IndelPod pod = pods.pods[p];
        if (pod.insert)
        {
            for (size_t i = 0; i < pod.size; i++)
            {
                size_t ins_pos = get_lateral(pod.indels[i].lateral3, insert3);
                if (ins_pos < ref_len)
                {
                    rels[ins_pos] |= ins_rel;
                }
            }
        }
    }
    
    // FIXME: AMBIGUOUS INDELS

    // Clip bases from the 5' and 3' ends of the read.
    read->ref_end5 = min(read->pos + clip_end5, ref_len + 1);
    read->ref_end3 = (ref_pos > clip_end3) ? (ref_pos - clip_end3) : 0;
    return 0;
}


static int calc_rels_line(unsigned char *rels,
                          SamRead *read,
                          const char *line,
                          const char *ref,
                          const char *ref_seq,
                          size_t ref_len,
                          unsigned long min_mapq,
                          unsigned char min_qual,
                          int insert3,
                          int ambindel,
                          size_t clip_end5,
                          size_t clip_end3)
{
    // Validate the SAM file line and parse it into a SamRead.
    assert(line != NULL);
    if (!(*line)) {return 1;}
    if (parse_sam_line(read, line)) {return 1;}
    if (validate_read(read, ref, min_mapq)) {return 1;}

    // Calculate relationships for the read.
    return calc_rels_read(rels,
                          read,
                          ref_seq,
                          ref_len,
                          min_qual,
                          insert3,
                          ambindel,
                          clip_end5,
                          clip_end3);
}


static int set_rel(PyObject *rels_dict, size_t pos, unsigned char rel)
{
    // Convert the C variables pos and rel into Python int objects.
    PyObject *key = PyLong_FromSize_t(pos);
    PyObject *value = PyLong_FromUnsignedLong((unsigned long)rel);

    if (key == NULL || value == NULL)
    {
        // Decrement reference counts if not NULL.
        Py_XDECREF(key);
        Py_XDECREF(value);
        Py_DECREF(rels_dict);
        return 1;
    }

    // Add the key-value pair to the dictionary.
    if (PyDict_SetItem(rels_dict, key, value) < 0)
    {
        Py_DECREF(key);
        Py_DECREF(value);
        Py_DECREF(rels_dict);
        return 1;
    }

    // Clean up the references for the key and value (they are now owned
    // by the dictionary).
    Py_DECREF(key);
    Py_DECREF(value);

    return 0;
}


static int put_rel_in_dict(PyObject *rels_dict,
                           size_t pos,
                           unsigned char rel)
{
    // Matches are not added to the rels_dict.
    if (rel != MATCH)
    {
        // It is impossible for any relationship within the region to be
        // non-covered (255, or 2^8 - 1).
        assert(rel != 255);
        if (set_rel(rels_dict, pos, rel))
            {return 1;}
    }
    return 0;
}


static int put_rels_in_dict(PyObject *rels_dict,
                            unsigned char *rels,
                            size_t end5,
                            size_t end3)
{
    // Validate the arguments; positions are 1-indexed and cannot be 0.
    assert(rels != NULL);
    if (end5 == 0 || end3 == 0)
    {
        PyErr_SetString(PyExc_ValueError, "Positions cannot be 0");
        return 1;
    }

    for (size_t pos = end5; pos <= end3; pos++)
    {
        // Subtract 1 from the 1-indexed position to make it a 0-indexed
        // array index.
        unsigned char rel = rels[pos - 1];
        // It is impossible for any relationship within the region to be
        // irreconcilable (0) for one read: only between two reads.
        assert(rel != 0);
        if (put_rel_in_dict(rels_dict, pos, rel))
            {return 1;}
    }
    return 0;
}


static int put_2_rels_in_dict(PyObject *rels_dict,
                              size_t fwd_end5,
                              size_t fwd_end3,
                              unsigned char *fwd_rels,
                              size_t rev_end5,
                              size_t rev_end3,
                              unsigned char *rev_rels)
{
    // Validate the arguments; positions are 1-indexed and cannot be 0.
    if (fwd_end5 == 0 || fwd_end3 == 0 || rev_end5 == 0 || rev_end3 == 0)
    {
        PyErr_SetString(PyExc_ValueError, "Positions cannot be 0");
        return 1;
    }

    // Find the region where both reads 1 and 2 overlap.
    // All positions are 1-indexed.
    size_t both_end5, both_end3;

    // Find the region 5' of the overlap.
    if (fwd_end5 > rev_end5)
    {
        // The forward read begins after the reverse read.
        both_end5 = fwd_end5;
        if (put_rels_in_dict(rels_dict,
                             rev_rels,
                             rev_end5,
                             min(rev_end3, both_end5 - 1)))
            {return 1;}
    }
    else
    {
        // The forward read begins with or after the reverse read.
        both_end5 = rev_end5;
        if (put_rels_in_dict(rels_dict,
                             fwd_rels,
                             fwd_end5,
                             min(fwd_end3, both_end5 - 1)))
            {return 1;}
    }

    // Find the region 3' of the overlap.
    if (fwd_end3 > rev_end3)
    {
        // The forward read ends after the reverse read.
        both_end3 = rev_end3;
        if (put_rels_in_dict(rels_dict,
                             fwd_rels,
                             fwd_end3,
                             max(fwd_end5, both_end3 + 1)))
            {return 1;}
    }
    else
    {
        // The forward read ends with or after the reverse read.
        both_end3 = fwd_end3;
        if (put_rels_in_dict(rels_dict,
                             rev_rels,
                             rev_end3,
                             max(rev_end5, both_end3 + 1)))
            {return 1;}
    }

    // Fill relationships in the region of overlap.
    for (size_t pos = both_end5; pos <= both_end3; pos++)
    {
        if (put_rel_in_dict(rels_dict, pos, fwd_rels[pos - 1] & rev_rels[pos - 1]))
            {return 1;}
    }

    return 0;
}


static int put_end_in_list(PyObject *ends_list, Py_ssize_t index, size_t end)
{
    PyObject *py_end = PyLong_FromSize_t(end);
    if (py_end == NULL) {return 1;}
    if (PyList_SetItem(ends_list, index, py_end)) {return 1;}
    return 0;
}


/*
Clean up by freeing all the dynamically allocated memory and optionally
decrementing the reference count for py_object, to prevent memory leaks.
*/
static PyObject *cleanup(unsigned char **rels1_ptr,
                         unsigned char **rels2_ptr,
                         PyObject *py_object)
{
    // All of the pointers to pointers must be non-NULL.
    assert(rels1_ptr != NULL);
    assert(rels2_ptr != NULL);

    // Free rels1 and point it to NULL.
    unsigned char *rels1 = *rels1_ptr;
    if (rels1 != NULL)
    {
        free(rels1);
        // Now that its memory has been freed, make rels1 point to NULL.
        *rels1_ptr = NULL;
    }

    // Free rels2 and point it to NULL.
    unsigned char *rels2 = *rels2_ptr;
    if (rels2 != NULL)
    {
        free(rels2);
        // Now that its memory has been freed, make rels2 point to NULL.
        *rels2_ptr = NULL;
    }

    // Decrement the reference count for py_object.
    if (py_object != NULL)
    {
        Py_DECREF(py_object);
    }

    // Return NULL so that py_calc_rels_lines can call this:
    // return error_cleanup(...);
    return NULL;
}


/* Python interface to function calc_rels_lines */
static PyObject *py_calc_rels_lines(PyObject *self, PyObject *args)
{
    // Attempt to convert the arguments from the Python function call
    // into C data types. Return NULL upon failure.
    const char *line1;
    const char *line2;
    const char *ref;
    const char *ref_seq;
    unsigned long ref_len;
    unsigned long min_mapq;
    unsigned char min_qual;
    int insert3;
    int ambindel;
    int overhangs;
    unsigned long clip_end5;
    unsigned long clip_end3;
    if (!PyArg_ParseTuple(args,
                          "sssskkbpppkk",
                          &line1,
                          &line2,
                          &ref,
                          &ref_seq,
                          &ref_len,
                          &min_mapq,
                          &min_qual,
                          &insert3,
                          &ambindel,
                          &overhangs,
                          &clip_end5,
                          &clip_end3))
        {return NULL;}
    
    // It is impossible for any of these pointers to be NULL.
    assert(line1 != NULL);
    assert(line2 != NULL);
    assert(ref != NULL);
    assert(ref_seq != NULL);

    // Determine the number of mates.
    Py_ssize_t num_mates = (*line2) ? 2 : 1;

    // Initialize containers to hold the results.
    unsigned char *rels1 = NULL, *rels2 = NULL;
    
    PyObject *ends_rels_tuple = PyTuple_New(2);
    if (ends_rels_tuple == NULL)
        {return cleanup(&rels1, &rels2, NULL);}
    
    PyObject *ends_tuple = PyTuple_New(2);
    if (ends_tuple == NULL)
        {return cleanup(&rels1, &rels2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_rels_tuple, 0, ends_tuple);

    PyObject *end5s_list = PyList_New(num_mates);
    if (end5s_list == NULL)
        {return cleanup(&rels1, &rels2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_tuple, 0, end5s_list);

    PyObject *end3s_list = PyList_New(num_mates);
    if (end3s_list == NULL)
        {return cleanup(&rels1, &rels2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_tuple, 1, end3s_list);
    
    PyObject *rels_dict = PyDict_New();
    if (rels_dict == NULL)
        {return cleanup(&rels1, &rels2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_rels_tuple, 1, rels_dict);

    // Allocate relationships for line 1.
    rels1 = calloc(ref_len, sizeof(rels1));
    if (rels1 == NULL)
    {
        PyErr_NoMemory();
        return cleanup(&rels1, &rels2, ends_rels_tuple);
    }

    // Calculate relationships for line 1.
    SamRead read1;
    if (calc_rels_line(rels1,
                       &read1,
                       line1,
                       ref,
                       ref_seq,
                       ref_len,
                       min_mapq,
                       min_qual,
                       insert3,
                       ambindel,
                       clip_end5,
                       clip_end3))
        {return cleanup(&rels1, &rels2, ends_rels_tuple);}

    if (num_mates > 1)
    {
        // The read comprises forward (fwd) and reverse (rev) mates.
        size_t fwd_end5, fwd_end3, rev_end5, rev_end3;
        // Check if line 2 differs from line 1.
        if (strcmp(line1, line2))
        {
            // Allocate relationships for line 2.
            rels2 = calloc(ref_len, sizeof(rels2));
            if (rels2 == NULL)
            {
                PyErr_NoMemory();
                return cleanup(&rels1, &rels2, ends_rels_tuple);
            }

            // Calculate relationships for line 2.
            SamRead read2;
            if (calc_rels_line(rels2,
                               &read2,
                               line2,
                               ref,
                               ref_seq,
                               ref_len,
                               min_mapq,
                               min_qual,
                               insert3,
                               ambindel,
                               clip_end5,
                               clip_end3))
                {return cleanup(&rels1, &rels2, ends_rels_tuple);}

            // Check if reads 1 and 2 are paired properly.
            if (validate_pair(&read1, &read2))
                {return cleanup(&rels1, &rels2, ends_rels_tuple);}

            // Determine the 5'/3' ends of the forward/reverse reads.
            unsigned char *fwd_rels = NULL, *rev_rels = NULL;
            if (read2.reverse)
            {
                fwd_end5 = read1.ref_end5;
                fwd_end3 = read1.ref_end3;
                fwd_rels = rels1;
                rev_end5 = read2.ref_end5;
                rev_end3 = read2.ref_end3;
                rev_rels = rels2;
            }
            else
            {
                fwd_end5 = read2.ref_end5;
                fwd_end3 = read2.ref_end3;
                fwd_rels = rels2;
                rev_end5 = read1.ref_end5;
                rev_end3 = read1.ref_end3;
                rev_rels = rels1;
            }

            // Remove overhangs if needed.
            if (!overhangs)
            {
                // The 5' end of the reverse mate cannot extend past
                // the 5' end of the forward mate.
                if (rev_end5 < fwd_end5) {rev_end5 = fwd_end5;}
                // The 3' end of the forward mate cannot extend past
                // the 3' end of the reverse mate.
                if (fwd_end3 > rev_end3) {fwd_end3 = rev_end3;}
            }

            // Merge reads 1 and 2.
            if (put_2_rels_in_dict(rels_dict,
                                   fwd_end5,
                                   fwd_end3,
                                   fwd_rels,
                                   rev_end5,
                                   rev_end3,
                                   rev_rels))
                {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        }
        else
        {
            // Lines 1 and 2 are identical.
            fwd_end5 = read1.ref_end5;
            fwd_end3 = read1.ref_end3;
            rev_end5 = read1.ref_end5;
            rev_end3 = read1.ref_end3;
            if (put_rels_in_dict(rels_dict,
                                 rels1,
                                 read1.ref_end5,
                                 read1.ref_end3))
                {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        }

        // Fill int the lists of 5' and 3' ends.
        if (put_end_in_list(end5s_list, 0, fwd_end5))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        if (put_end_in_list(end5s_list, 1, rev_end5))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 0, fwd_end3))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 1, rev_end3))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
    }
    else
    {
        // Line 2 does not exist.
        if (put_rels_in_dict(rels_dict,
                             rels1,
                             read1.ref_end5,
                             read1.ref_end3))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        
        // Fill int the lists of 5' and 3' ends.
        if (put_end_in_list(end5s_list, 0, read1.ref_end5))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 0, read1.ref_end3))
            {return cleanup(&rels1, &rels2, ends_rels_tuple);}
    }

    // Free the memory of rels1 and rels2, but not of the object that
    // will be returned.
    cleanup(&rels1, &rels2, NULL);

    return ends_rels_tuple;
}


/* Python method table */
static PyMethodDef RelateMethods[] = {
    {"calc_rels_lines",
     py_calc_rels_lines,
     METH_VARARGS,
     "Calculate the relationships for a line or pair of lines"},
    {NULL, NULL, 0, NULL}  // sentinel
};


/* Python module definition */
static struct PyModuleDef relatemodule = {
    PyModuleDef_HEAD_INIT,
    "relate",      // module name
    NULL,          // TODO documentation
    -1,            // module state
    RelateMethods  // method table
};


/* Python module initialization function */
PyMODINIT_FUNC PyInit_relate(void)
{
    // Initialize the module and ensure it exists.
    PyObject *module;
    module = PyModule_Create(&relatemodule);
    if (module == NULL) {return NULL;}

    // Define a new type of Python exception.
    RelateError = PyErr_NewException("relate.RelateError", NULL, NULL);
    // Add the exception type to the module.
    Py_XINCREF(RelateError);
    if (PyModule_AddObject(module, "RelateError", RelateError) < 0)
    {
        // Adding the exception type failed. Stop and return NULL.
        Py_XDECREF(RelateError);
        Py_CLEAR(RelateError);
        Py_DECREF(module);
        return NULL;
    }

    // Initializing the module succeeded.
    return module;
}


int main()
{
    // This module must be called via its Python API, not run directly.
    return 0;
}
