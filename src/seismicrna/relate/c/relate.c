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
    0 if successful, otherwise -1
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
        return -1;
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
    0 if successful, otherwise -1
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
        return -1;
    }

    // Bitwise flag
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end)))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse SAM flag");
        return -1;
    }
    read->flag = temp_ulong;
    if (read->flag > MAX_FLAG)
    {
        PyErr_SetString(PyExc_ValueError, "SAM flag is too large");
        return -1;
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
        return -1;
    }

    // Mapping position
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end)))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse mapping position");
        return -1;
    }
    read->pos = (size_t)temp_ulong;

    // Mapping quality
    if (parse_ulong(&temp_ulong, strtok_r(NULL, SAM_SEP, &end)))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse mapping quality");
        return -1;
    }
    read->mapq = temp_ulong;

    // CIGAR string
    if ((read->cigar = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse CIGAR string");
        return -1;
    }

    // Next reference (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse next reference name");
        return -1;
    }

    // Next position (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse next mapping position");
        return -1;
    }

    // Template length (ignored)
    if (strtok_r(NULL, SAM_SEP, &end) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse template length");
        return -1;
    }

    // Read sequence
    if ((read->seq = strtok_r(NULL, SAM_SEP, &end)) == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse read sequence");
        return -1;
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
        return -1;
    }

    // Lengths of read and quality strings must match.
    // Subtract 1 from end for the same reason as in "Read end position"
    if (end != NULL && read->len != (size_t)((end - 1) - read->qual))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Sequence and quality strings differ in length");
        return -1;
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
        return -1;
    }
    if (strcmp(read->ref, ref))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Reference name does not match name of SAM file");
        return -1;
    }
    return 0;
}


static int validate_pair(const SamRead *read1,
                         const SamRead *read2)
{
    // None of the pointers can be NULL.
    assert(read1 != NULL);
    assert(read2 != NULL);
    assert(read1->name != NULL);
    assert(read2->name != NULL);

    // Mates 1 and 2 must have the same name.
    if (strcmp(read1->name, read2->name))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 have different names");
        return -1;
    }

    // Both mates must be paired-end.
    if (!(read1->paired && read2->paired))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 are not both paired-end");
        return -1;
    }

    // Mates 1 and 2 must be marked as first and second, respectively.
    if (!(read1->first && read2->second))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 are not marked as first and second");
        return -1;
    }

    // Mates 1 and 2 must be in opposite orientations (i.e. one must be
    // marked as reversed and the other as not reversed).
    if (read1->reverse == read2->reverse)
    {
        PyErr_SetString(PyExc_ValueError,
                        "Mates 1 and 2 aligned in the same orientation");
        return -1;
    }

    return 0;
}


// CIGAR string operations

typedef struct
{
    const char *op;  // type of operation
    size_t len;      // length of operation
} CigarOp;


/*
Find the next operation in a CIGAR string. Store the kind of operation
in cigar->op and the length of the operation in cigar->len. If parsing
does not return a valid operation (which must have cigar->len > 0 and
cigar->op != NULL), then set cigar->op to NULL. This event happens when
the end of the CIGAR string is reached and does not signify an error.

Parameters
----------
cigar
    Non-nullable pointer to a CIGAR operation struct in which to store
    the parsed information from the CIGAR string

Returns
-------
int
    0 if successful, otherwise -1
*/
static inline int get_next_cigar_op(CigarOp *cigar)
{
    // Parse as many characters of text as possible to a number, cast
    // to size_t, and store in cigar->len (length of CIGAR operation).
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
        if (*cigar->op)
        {
            PyErr_SetString(PyExc_ValueError,
                            "Failed to parse a CIGAR operation");
            return -1;
        }
        // Point cigar->op to NULL to signify reaching the end.
        cigar->op = NULL;
    }
    return 0;
}


/* Encode a base character as a substitution. */
static inline unsigned char encode_subs(char base)
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


/* Encode a match or substitution as such if the quality is sufficient,
otherwise as ambiguous. */
static inline unsigned char encode_relate(char ref_base,
                                          char read_base,
                                          char read_qual,
                                          unsigned char min_qual)
{
    if (read_qual < min_qual) {return ANY_N ^ encode_subs(ref_base);}
    return (ref_base == read_base) ? MATCH : encode_subs(read_base);
}


/* More efficient version of encode_relate that assumes the base call
is not a substitution. */
static inline unsigned char encode_match(char read_base,
                                         char read_qual,
                                         unsigned char min_qual)
{
    return (read_qual >= min_qual) ? MATCH : (ANY_N ^ encode_subs(read_base));
}


/* Return whether two relationships are consistent with each other. */
static inline int consistent_rels(char rel1, char rel2)
{
    // Two relationships are consistent if they have one or more primary
    // relationships (bits) in common.
    if (rel1 & rel2) {return 1;}
    // They are also consistent if they are both substitutions.
    if ((rel1 & SUB_N) && (rel2 & SUB_N)) {return 1;}
    // Otherwise, they are inconsistent.
    return 0;
}


// Insertions and deletions

static const int DELETION = 0;
static const int INSERTION = 1;
static const int INIT_PODS_CAPACITY = 4;
static const size_t INIT_POD_CAPACITY = 4;
static const size_t CAPACITY_FACTOR = 2;


static inline unsigned char get_ins_rel(int insert3)
{
    return insert3 ? INS_3 : INS_5;
}


static inline size_t calc_lateral5(size_t lateral3)
{
    return (lateral3 > 0) ? lateral3 - 1 : 0;
}


static inline size_t get_lateral(size_t lateral3, int insert3)
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
    size_t order;
    int insert;
    Indel *indels;
    size_t capacity;
    size_t num_indels;
} IndelPod;


typedef struct
{
    IndelPod *pods;
    size_t capacity;
    size_t num_pods;
} IndelPodArray;


static int init_pod(IndelPod *pod, size_t order, int insert)
{
    pod->order = order;
    pod->insert = insert;
    // Assume that the pod is being initialized because an indel needs
    // to go into the pod, so initialize with the capacity to hold up to
    // INIT_POD_CAPACITY indels.
    pod->capacity = INIT_POD_CAPACITY;
    pod->indels = malloc(pod->capacity * sizeof(Indel));
    if (pod->indels == NULL)
    {
        PyErr_NoMemory();
        return -1;
    }
    // printf("malloc %p\n", pod->indels);
    pod->num_indels = 0;
    return 0;
}


static int add_indel(IndelPodArray *pods,
                     int insert,
                     size_t opposite,
                     size_t lateral3)
{
    IndelPod *pod = NULL;
    Indel *indel = NULL;
    if (pods->num_pods == 0)
    {
        assert(pods->pods == NULL);
        assert(pods->capacity == 0);
        // There are no pods yet: allocate a new array of pods with the
        // capacity to hold up to INIT_PODS_CAPACITY pods.
        pods->capacity = INIT_PODS_CAPACITY;
        pods->pods = malloc(pods->capacity * sizeof(IndelPod));
        if (pods->pods == NULL)
        {
            PyErr_NoMemory();
            return -1;
        }
        // printf("malloc %p\n", pods->pods);
        // Initialize the first pod in the new array.
        pod = pods->pods;
        if (init_pod(pod, 0, insert)) {return -1;}
        // Increment num_pods only if init_pod succeeds because, if it
        // fails, then pod->indels is NULL, so free_pods() would try to
        // free() a NULL pointer if num_pods were already incremented.
        pods->num_pods = 1;
    }
    else
    {
        assert(pods->pods != NULL);
        assert(pods->capacity >= pods->num_pods);
        // There are pods: check if the last pod is the correct type.
        pod = &(pods->pods[pods->num_pods - 1]);
        if (pod->insert != insert)
        {
            // The last pod is the other type: a new pod is needed.
            if (pods->num_pods == pods->capacity)
            {
                // The array of pods is already at its maximum capacity:
                // allocate a new array with twice the capacity.
                pods->capacity *= CAPACITY_FACTOR;
                pod = realloc(pods->pods, pods->capacity * sizeof(IndelPod));
                if (pod == NULL)
                {
                    PyErr_NoMemory();
                    return -1;
                }
                // printf("implicit-free %p\n", pods->pods);
                // printf("realloc %p\n", pod);
                pods->pods = pod;
            }
            // Initialize the new pod.
            pod = &(pods->pods[pods->num_pods]);
            if (init_pod(pod, pods->num_pods, insert)) {return -1;}
            // Increment num_pods only if init_pod succeeds because, if
            // it fails, then pod->indels is NULL, so free_pods() would
            // try to free() a NULL pointer if num_pods were already
            // incremented.
            pods->num_pods++;
        }
        else
        {
            // The last pod is the correct type.
            if (pod->num_indels == pod->capacity)
            {
                // The pod is already at its maximum capacity: allocate
                // a new array with twice the capacity.
                pod->capacity *= CAPACITY_FACTOR;
                indel = realloc(pod->indels, pod->capacity * sizeof(Indel));
                if (indel == NULL)
                {
                    PyErr_NoMemory();
                    return -1;
                }
                // printf("implicit-free %p\n", pod->indels);
                // printf("realloc %p\n", indel);
                pod->indels = indel;
            }
        }
    }
    // Initialize the new indel.
    assert(pod != NULL);
    assert(pod->num_indels < pod->capacity);
    indel = &(pod->indels[pod->num_indels]);
    indel->insert = insert;
    indel->opposite = opposite;
    indel->lateral3 = lateral3;
    pod->num_indels++;
    return 0;
}


static void free_pods(IndelPodArray *pods)
{
    assert(pods != NULL);
    // Free the indel pods.
    if (pods->num_pods > 0)
    {
        assert(pods->pods != NULL);
        // Free the indels in each pod.
        for (size_t p = 0; p < pods->num_pods; p++)
        {
            assert(pods->pods[p].indels != NULL);
            // printf("free %p\n", pods->pods[p].indels);
            free(pods->pods[p].indels);
        }
        // printf("free %p\n", pods->pods);
        free(pods->pods);
        // Now that its memory has been freed, point pods->pods at NULL.
        pods->pods = NULL;
        pods->num_pods = 0;
    }
    else
    {
        assert(pods->pods == NULL);
    }
}


/* Determine the order of two IndelPods. */
static int comp_pods_fwd(const void *a, const void *b)
{
    const IndelPod *pod1 = (const IndelPod *)a;
    const IndelPod *pod2 = (const IndelPod *)b;
    // Two IndelPods cannot have the same value for order.
    assert(pod1->order != pod2->order);
    // Use comparison, not subtraction (i.e. pod1->order - pod2->order),
    // because order is size_t (which is unsigned); if the difference is
    // negative, then it will overflow to a very large positive number.
    return pod1->order > pod2->order ? 1 : -1;
}


/* Determine the reverse order of two IndelPods. */
static int comp_pods_rev(const void *a, const void *b)
{
    const IndelPod *pod1 = (const IndelPod *)a;
    const IndelPod *pod2 = (const IndelPod *)b;
    // Two IndelPods cannot have the same value for order.
    assert(pod1->order != pod2->order);
    // Use comparison, not subtraction (i.e. pod1->order - pod2->order),
    // because order is size_t (which is unsigned); if the difference is
    // negative, then it will overflow to a very large positive number.
    return pod1->order < pod2->order ? 1 : -1;
}


/* Sort the IndelPods in an IndelPodArray, either forward or reverse. */
static void sort_pods(IndelPodArray *pods, int move5to3)
{
    if (move5to3)
        {qsort(pods->pods, pods->num_pods, sizeof(IndelPod), comp_pods_fwd);}
    else
        {qsort(pods->pods, pods->num_pods, sizeof(IndelPod), comp_pods_rev);}
}


/* Determine the order of two Indels. */
static int comp_indels_fwd(const void *a, const void *b)
{
    const Indel *indel1 = (const Indel *)a;
    const Indel *indel2 = (const Indel *)b;
    // Two Indels cannot have the same value for opposite.
    assert(indel1->opposite != indel2->opposite);
    // Use comparison, not subtraction, because opposite is size_t
    // (which is unsigned); if the difference is negative, then it will
    // overflow to a very large positive number.
    return indel1->opposite > indel2->opposite ? 1 : -1;
}


/* Determine the reverse order of two Indels. */
static int comp_indels_rev(const void *a, const void *b)
{
    const Indel *indel1 = (const Indel *)a;
    const Indel *indel2 = (const Indel *)b;
    // Two Indels cannot have the same value for opposite.
    assert(indel1->opposite != indel2->opposite);
    // Use comparison, not subtraction, because opposite is size_t
    // (which is unsigned); if the difference is negative, then it will
    // overflow to a very large positive number.
    return indel1->opposite < indel2->opposite ? 1 : -1;
}


/* Sort the Indels in an IndelPod, either forward or reverse. */
static void sort_pod(IndelPod *pod, int move5to3)
{
    if (move5to3)
        {qsort(pod->indels, pod->num_indels, sizeof(Indel), comp_indels_fwd);}
    else
        {qsort(pod->indels, pod->num_indels, sizeof(Indel), comp_indels_rev);}
}


/* Move one indel. */
static inline void move_indel(Indel *indel, size_t opposite, size_t lateral3)
{
    indel->opposite = opposite;
    indel->lateral3 = lateral3;
}


/* Move one indel while adjusting the positions of any other indels in
its pod through which the moving indel tunnels. */
static void move_indels(IndelPod *pod,
                        size_t indel_index,
                        size_t opposite,
                        size_t lateral3)
{
    // Determine which indel is at the given index.
    assert(indel_index < pod->num_indels);
    Indel *indel = &(pod->indels[indel_index]);
    // Move every other indel in the pod that lies between the given
    // indel and the position to which the given indel should be moved.
    for (size_t i = 0; i < pod->num_indels; i++)
    {
        if (i != indel_index)
        {
            // Check if the other indel lies between the given indel and
            // the position to which the given indel should be moved.
            Indel *other = &(pod->indels[i]);
            if ((indel->opposite < other->opposite &&
                 other->opposite < opposite) ||
                (indel->opposite > other->opposite &&
                 other->opposite > opposite))
            {
                // If so, then move it.
                move_indel(other, opposite, lateral3);
            }
        }
    }
    // Move the given indel.
    move_indel(indel, opposite, lateral3);
}


/* Return whether a pod has an indel with the given opposite value. */
static inline int pod_has_opposite(IndelPod *pod, size_t opposite)
{
    for (size_t i = 0; i < pod->num_indels; i++)
    {
        if (pod->indels[i].opposite == opposite) {return 1;}
    }
    return 0;
}


/* Calculate the positions for moving an indel. */
static void calc_positions(size_t *swap_lateral,
                          size_t *next_lateral3,
                          size_t *next_opposite,
                          IndelPod *pod,
                          Indel *indel,
                          int move5to3)
{
    if (move5to3)
    {
        // Find the position within the indel's own sequence with which
        // to try to swap the indel.
        *swap_lateral = indel->lateral3;
        // If the indel moves, then its position within its own sequence
        // will change.
        *next_lateral3 = indel->lateral3 + 1;
        // Find the position within the opposite sequence to which to
        // try to move the indel.
        *next_opposite = indel->opposite + 1;
        while (pod_has_opposite(pod, *next_opposite))
            {(*next_opposite)++;}
    }
    else
    {
        // Find the position within the indel's own sequence with which
        // to try to swap the indel.
        *swap_lateral = calc_lateral5(indel->lateral3);
        // If the indel moves, then its position within its own sequence
        // will change.
        assert(indel->lateral3 > 0);
        *next_lateral3 = indel->lateral3 - 1;
        // Find the position within the opposite sequence to which to
        // try to move the indel.
        assert(indel->opposite > 0);
        *next_opposite = indel->opposite - 1;
        while (pod_has_opposite(pod, *next_opposite))
            assert(*next_opposite > 0);
            {(*next_opposite)--;}
    }
}


/* Check whether position is out of bounds. */
static inline int check_out_of_bounds(size_t position,
                                      size_t ref_end5,
                                      size_t ref_end3)
{
    // 
}


/* Try to move a deletion. */
static int try_move_del(unsigned char *rels,
                        IndelPodArray *pods,
                        SamRead *read,
                        const char *ref_seq,
                        size_t ref_len,
                        size_t ref_end5,
                        size_t ref_end3,
                        size_t read_end5,
                        size_t read_end3,
                        unsigned char min_qual,
                        int insert3,
                        int move5to3,
                        size_t pod_index,
                        size_t indel_index)
{
    assert(pod_index < pods->num_pods);
    IndelPod *pod = &(pods->pods[pod_index]);
    assert(indel_index < pod->num_indels);
    Indel *indel = &(pod->indels[indel_index]);

    // Calculate the positions to which to move the deletion.
    size_t swap_lateral;
    size_t next_lateral3;
    size_t next_opposite;
    calc_positions(&swap_lateral,
                   &next_lateral3,
                   &next_opposite,
                   pod,
                   indel,
                   move5to3);
}


static int try_move_indel(unsigned char *rels,
                          IndelPodArray *pods,
                          SamRead *read,
                          const char *ref_seq,
                          size_t ref_len,
                          size_t ref_end5,
                          size_t ref_end3,
                          size_t read_end5,
                          size_t read_end3,
                          unsigned char min_qual,
                          int insert3,
                          int move5to3,
                          size_t pod_index,
                          size_t indel_index)
{

}


static void find_ambindels_recurse(unsigned char *rels,
                                   IndelPodArray *pods,
                                   SamRead *read,
                                   const char *ref_seq,
                                   size_t ref_len,
                                   size_t ref_end5,
                                   size_t ref_end3,
                                   size_t read_end5,
                                   size_t read_end3,
                                   unsigned char min_qual,
                                   int insert3,
                                   int move5to3,
                                   size_t pod_index,
                                   size_t indel_index)
{
    // Sort the indels in the pod and select the indel at this index.
    assert(pod_index < pods->num_pods);
    IndelPod *pod = &(pods->pods[pod_index]);
    sort_pod(pod, move5to3);
    assert(indel_index < pod->num_indels);
    Indel *indel = &(pod->indels[indel_index]);

    /*
    // Record the initial position of the indel at this index.
    size_t init_opposite = indel->opposite;
    size_t init_lateral3 = indel->lateral3;

    // Try to move the indel at this index one step.
    if (try_move_indel())
    {
        // Try to move the indel at this index one more step.
        find_ambindels_recurse(rels,
                               pods,
                               read,
                               ref_seq,
                               ref_len,
                               ref_end5,
                               ref_end3,
                               read_end5,
                               read_end3,
                               min_qual,
                               insert3,
                               move5to3,
                               pod_index,
                               indel_index);
        if (move5to3)
        {
            // Move the indel back to its initial position.
            move_indels(pod, indel_index, init_opposite, init_lateral3);
        }
    }
    */

    // Check if the pod contains another indel after the current one.
    if (indel_index + 1 < pod->num_indels)
    {
        // Move the next indel in the pod.
        find_ambindels_recurse(rels,
                               pods,
                               read,
                               ref_seq,
                               ref_len,
                               ref_end5,
                               ref_end3,
                               read_end5,
                               read_end3,
                               min_qual,
                               insert3,
                               move5to3,
                               pod_index,
                               indel_index + 1);
    }

    // Check if there is another pod before the current one.
    if (pod_index > 0)
    {
        // Move the indels in the previous pod.
        find_ambindels_recurse(rels,
                               pods,
                               read,
                               ref_seq,
                               ref_len,
                               ref_end5,
                               ref_end3,
                               read_end5,
                               read_end3,
                               min_qual,
                               insert3,
                               move5to3,
                               pod_index - 1,
                               0);
    }
}


static void find_ambindels(unsigned char *rels,
                           IndelPodArray *pods,
                           SamRead *read,
                           const char *ref_seq,
                           size_t ref_len,
                           size_t ref_end5,
                           size_t ref_end3,
                           size_t read_end5,
                           size_t read_end3,
                           unsigned char min_qual,
                           int insert3)
{
    if (pods->num_pods > 0)
    {
        assert(pods->pods != NULL);
        for (int move5to3 = 0; move5to3 <= 1; move5to3++)
        {
            sort_pods(pods, move5to3);
            find_ambindels_recurse(rels,
                                   pods,
                                   read,
                                   ref_seq,
                                   ref_len,
                                   ref_end5,
                                   ref_end3,
                                   read_end5,
                                   read_end3,
                                   min_qual,
                                   insert3,
                                   move5to3,
                                   pods->num_pods - 1,
                                   0);
        }
    }
}


/* Count the reference bases consumed by the CIGAR string. */
static int count_ref_bases(const char *read_cigar, size_t *ref_bases)
{
    assert(read_cigar != NULL);
    assert(ref_bases != NULL);

    // Initialize the count to 0.
    *ref_bases = 0;

    // Initialize the CIGAR operation so that it points just before the
    // CIGAR string.
    CigarOp cigar;
    cigar.op = read_cigar - 1;

    // Read the first operation from the CIGAR string; catch errors.
    if (get_next_cigar_op(&cigar)) {return -1;}

    // Read the entire CIGAR string one operation at a time.
    while (cigar.op != NULL)
    {
        // Decide what to do based on the current CIGAR operation.
        switch (*cigar.op)
        {
            case CIG_ALIGN:
            case CIG_MATCH:
            case CIG_SUBST:
            case CIG_DELET:
                // These operations consume the reference.
                *ref_bases += cigar.len;
                break;
            default:
                // Other operations do not.
                break;
        }

        // Read the next CIGAR operation; catch errors.
        if (get_next_cigar_op(&cigar)) {return -1;}
    }

    return 0;
}


/* Compute the relationships of a SamRead. */
static int calc_rels_read(unsigned char *rels,
                          IndelPodArray *pods,
                          const SamRead *read,
                          const char *rels_seq,
                          size_t rels_len,
                          unsigned char min_qual,
                          int insert3,
                          int ambindel)
{
    // Validate the arguments.
    assert(rels != NULL);
    assert(pods != NULL);
    assert(read != NULL);
    assert(rels_seq != NULL);

    // Positions in the relationship array and read (0-indexed).
    size_t rels_pos = 0;
    size_t read_pos = 0;
    // 5' and 3' ends of the read, in the read coordinates (1-indexed).
    size_t read_end5 = 1;
    size_t read_end3 = read->len;
    // CIGAR operation stopping point.
    size_t cigar_op_stop_pos;
    
    // Initialize the CIGAR operation so that it points just before the
    // CIGAR string.
    assert(read->cigar != NULL);
    CigarOp cigar;
    cigar.op = read->cigar - 1;

    // Read the first operation from the CIGAR string; catch errors.
    if (get_next_cigar_op(&cigar)) {return -1;}

    // Read the entire CIGAR string one operation at a time.
    while (cigar.op != NULL)
    {
        // Decide what to do based on the current CIGAR operation.
        switch (*cigar.op)
        {
    
            case CIG_MATCH:
                // The read and reference match.
                cigar_op_stop_pos = rels_pos + cigar.len;
                if (cigar_op_stop_pos > rels_len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A match extended out of the reference"
                    );
                    return -1;
                }
                if (read_pos + cigar.len > read->len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A match extended out of the read"
                    );
                    return -1;
                }
                while (rels_pos < cigar_op_stop_pos)
                {
                    rels[rels_pos] = encode_match(read->seq[read_pos],
                                                  read->qual[read_pos],
                                                  min_qual);
                    rels_pos++;
                    read_pos++;
                }
                break;
            
            case CIG_ALIGN:
            case CIG_SUBST:
                // The read and reference have matches or substitutions.
                cigar_op_stop_pos = rels_pos + cigar.len;
                if (cigar_op_stop_pos > rels_len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "An operation extended out of the reference"
                    );
                    return -1;
                }
                if (read_pos + cigar.len > read->len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "An operation extended out of the read"
                    );
                    return -1;
                }
                while (rels_pos < cigar_op_stop_pos)
                {
                    rels[rels_pos] = encode_relate(rels_seq[rels_pos],
                                                   read->seq[read_pos],
                                                   read->qual[read_pos],
                                                   min_qual);
                    rels_pos++;
                    read_pos++;
                }
                break;
            
            case CIG_DELET:
                // A portion of the reference is deleted from the read.
                cigar_op_stop_pos = rels_pos + cigar.len;
                if (cigar_op_stop_pos > rels_len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A deletion extended out of the reference"
                    );
                    return -1;
                }
                if (rels_pos == 0)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A deletion occured at the beginning of the reference"
                    );
                    return -1;
                }
                if (read_pos == 0)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A deletion occured at the beginning of the read"
                    );
                    return -1;
                }
                if (read_pos >= read->len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A deletion occured at the end of the read"
                    );
                    return -1;
                }
                while (rels_pos < cigar_op_stop_pos)
                {
                    // Mark the deletion in the array of relationships.
                    rels[rels_pos] = DELET;
                    // Add a deletion to the record of indels.
                    if (add_indel(pods, DELETION, rels_pos, read_pos))
                        {return -1;}
                    rels_pos++;
                }
                break;
            
            case CIG_INSRT:
                // The read contains an insertion.
                cigar_op_stop_pos = read_pos + cigar.len;
                if (cigar_op_stop_pos > read->len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "An insertion extended out of the read"
                    );
                    return -1;
                }
                if (read_pos == 0)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "An insertion occured at the beginning of the read"
                    );
                    return -1;
                }
                if (rels_pos == 0)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "An insertion occured at the beginning of the reference"
                    );
                    return -1;
                }
                if (rels_pos >= rels_len)
                {
                    // printf("REF POS: %zu\n", ref_pos);
                    // printf("REF LEN: %zu\n", ref_len);
                    PyErr_SetString(
                        PyExc_ValueError,
                        "An insertion occured at the end of the reference"
                    );
                    return -1;
                }
                while (read_pos < cigar_op_stop_pos)
                {
                    // Add an insertion to the record of indels.
                    if (add_indel(pods, INSERTION, read_pos, rels_pos))
                        {return -1;}
                    read_pos++;
                }
                break;
            
            case CIG_SCLIP:
                // Bases were soft-clipped from the 5' or 3' end of the
                // read during alignment. Like insertions, they consume
                // the read but not the reference; however, they are not
                // mutations.
                cigar_op_stop_pos = read_pos + cigar.len;
                if (cigar_op_stop_pos > read->len)
                {
                    PyErr_SetString(
                        PyExc_ValueError,
                        "A soft clip extended out of the read"
                    );
                    return -1;
                }
                if (read_pos == 0)
                {
                    // This is the 5' soft clip.
                    if (read_end5 != 1)
                    {
                        PyErr_SetString(
                            PyExc_ValueError,
                            "The read contained > 1 soft clip at the beginning"
                        );
                        return -1;
                    }
                    // Clip bases from the the 5' end of the read.
                    read_end5 += cigar.len;
                }
                else
                {
                    // This is the 3' soft clip.
                    if (cigar_op_stop_pos != read->len)
                    {
                        PyErr_SetString(
                            PyExc_ValueError,
                            "The read contained a soft clip in the middle"
                        );
                        return -1;
                    }
                    if (read_end3 != read->len)
                    {
                        PyErr_SetString(
                            PyExc_ValueError,
                            "The read contained > 1 soft clip at the end"
                        );
                        return -1;
                    }
                    // Clip bases from the 3' end of the read.
                    read_end3 -= cigar.len;
                }
                read_pos = cigar_op_stop_pos;
                break;
            
            default:
                // The CIGAR operation was not recognized.
                PyErr_SetString(
                    PyExc_ValueError,
                    "The CIGAR string contained an unknown type of operation"
                );
                return -1;
        }

        // Read the next operation from the CIGAR string; catch errors.
        if (get_next_cigar_op(&cigar)) {return -1;}
    }

    // Verify that the sum of the lengths of all CIGAR operations that
    // consumed the reference equals the length of the relationships.
    if (rels_pos != rels_len)
    {
        PyErr_SetString(PyExc_ValueError,
                        "The number of reference bases consumed by the CIGAR "
                        "string differed from the length of the relationships");
        return -1;
    }
    // Verify that the sum of the lengths of all CIGAR operations that
    // consumed the read equals the length of the read.
    if (read_pos != read->len)
    {
        PyErr_SetString(PyExc_ValueError,
                        "The number of read bases consumed by the CIGAR "
                        "string differed from the length of the read");
        return -1;
    }

    // Add insertions to rels.
    const unsigned char ins_rel = get_ins_rel(insert3);
    for (size_t p = 0; p < pods->num_pods; p++)
    {
        IndelPod *pod = &(pods->pods[p]);
        if (pod->insert)
        {
            for (size_t i = 0; i < pod->num_indels; i++)
            {
                size_t ins_pos = get_lateral(pod->indels[i].lateral3, insert3);
                if (ins_pos < rels_len)
                {
                    rels[ins_pos] |= ins_rel;
                }
            }
        }
    }
    
    if (ambindel && pods->num_pods > 0)
    {
        // Find and mark ambiguous insertions and deletions.
        find_ambindels(rels,
                       pods,
                       read,
                       rels_seq,
                       rels_len,
                       read->pos,
                       rels_pos,
                       read_end5,
                       read_end3,
                       min_qual,
                       insert3);
    }

    return 0;
}


static int calc_rels_line(unsigned char **rels,
                          IndelPodArray *pods,
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
    if (!(*line))
    {
        PyErr_SetString(PyExc_ValueError,
                        "Got an empty line to parse as a SAM read");
        return -1;
    }
    if (parse_sam_line(read, line)) {return -1;}
    if (validate_read(read, ref, min_mapq)) {return -1;}

    // Validate the read's start position.
    if (read->pos == 0)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "Mapping position is 0"
        );
        return -1;
    }
    if (read->pos > ref_len)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "Mapping position is greater than length of reference sequence"
        );
        return -1;
    }

    // Count the reference bases consumed by the CIGAR string.
    size_t ref_bases;
    if (count_ref_bases(read->cigar, &ref_bases))
        {return 1;}
    if (ref_bases == 0)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "The CIGAR string consumes 0 bases in the reference"
        );
        return -1;
    }

    // Calculate the reference position at which the 5' end of the read
    // is located, which is the mapping position plus the 5' clip (up to
    // a maximum of the reference length plus 1).
    read->ref_end5 = min(read->pos + clip_end5, ref_len + 1);
    // Calculate the reference position at which the 3' end of the read
    // is located, which is the mapping position plus one less than the
    // number of reference bases consumed, minus the 3' clip (down to a
    // minimum of 0).
    // Because ref_bases is guaranteed to be > 0, subtracting 1 will not
    // yield a negative number and cause overflow.
    read->ref_end3 = read->pos + (ref_bases - 1);
    // Also to avoid overflow, subtract clip_end3 only after confirming
    // that the difference would not be negative.
    read->ref_end3 =
        (read->ref_end3 > clip_end3) ? (read->ref_end3 - clip_end3) : 0;

    // Allocate memory for the relationships. The relationships will be
    // filled using the CIGAR string, so the number of relationships
    // allocated must equal the number of reference bases consumed by
    // the CIGAR string, not the number after clipping.
    assert(rels != NULL);
    assert(*rels == NULL);  // If *rels != NULL, then memory will leak.
    *rels = malloc(ref_bases * sizeof(*rels));
    if (*rels == NULL)
    {
        PyErr_NoMemory();
        return -1;
    }
    // printf("malloc %p\n", rels1);

    // Calculate relationships for the read.
    // Because read->pos is guaranteed to be > 0, subtracting 1 will not
    // yield a negative number and cause overflow.
    return calc_rels_read(*rels,
                          pods,
                          read,
                          ref_seq + (read->pos - 1),
                          ref_bases,
                          min_qual,
                          insert3,
                          ambindel);
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
        return -1;
    }

    // Add the key-value pair to the dictionary.
    if (PyDict_SetItem(rels_dict, key, value) < 0)
    {
        Py_DECREF(key);
        Py_DECREF(value);
        Py_DECREF(rels_dict);
        return -1;
    }

    // Clean up the references for the key and value (they are now owned
    // by the dictionary).
    Py_DECREF(key);
    Py_DECREF(value);

    return 0;
}


static inline int put_rel_in_dict(PyObject *rels_dict,
                                  size_t pos,
                                  unsigned char rel)
{
    // Matches are not added to the rels_dict.
    if (rel != MATCH)
    {
        // It is impossible for any relationship within the region to be
        // non-covered (255, or 2^8 - 1).
        assert(rel != 255);
        if (set_rel(rels_dict, pos, rel)) {return -1;}
    }
    return 0;
}


static int put_rels_in_dict(PyObject *rels_dict,
                            unsigned char *rels,
                            size_t read_pos,
                            size_t end5,
                            size_t end3)
{
    // Validate the arguments; positions are 1-indexed, and only end3
    // is allowed to be 0.
    assert(rels != NULL);
    assert(read_pos > 0);
    assert(end5 >= read_pos);

    for (size_t pos = end5; pos <= end3; pos++)
    {
        // Since rels is an array whose 0th index is the read's mapping
        // position, subtract the mapping position from the position of
        // the relationship to find its 0-based index in rels.
        unsigned char rel = rels[pos - read_pos];
        // It is impossible for any relationship within the region to be
        // irreconcilable (0) for one read: only between two reads.
        assert(rel != 0);
        if (put_rel_in_dict(rels_dict, pos, rel)) {return -1;}
    }
    return 0;
}


static int put_2_rels_in_dict(PyObject *rels_dict,
                              size_t fwd_read_pos,
                              size_t fwd_end5,
                              size_t fwd_end3,
                              unsigned char *fwd_rels,
                              size_t rev_read_pos,
                              size_t rev_end5,
                              size_t rev_end3,
                              unsigned char *rev_rels)
{
    // Validate the arguments; positions are 1-indexed and the 5' ends
    // must be > 0.
    assert(fwd_end5 > 0);
    assert(rev_end5 > 0);

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
                             rev_read_pos,
                             rev_end5,
                             min(rev_end3, both_end5 - 1)))
            {return -1;}
    }
    else
    {
        // The forward read begins with or after the reverse read.
        both_end5 = rev_end5;
        if (put_rels_in_dict(rels_dict,
                             fwd_rels,
                             fwd_read_pos,
                             fwd_end5,
                             min(fwd_end3, both_end5 - 1)))
            {return -1;}
    }

    // Find the region 3' of the overlap.
    if (fwd_end3 > rev_end3)
    {
        // The forward read ends after the reverse read.
        both_end3 = rev_end3;
        if (put_rels_in_dict(rels_dict,
                             fwd_rels,
                             fwd_read_pos,
                             fwd_end3,
                             max(fwd_end5, both_end3 + 1)))
            {return -1;}
    }
    else
    {
        // The forward read ends with or after the reverse read.
        both_end3 = fwd_end3;
        if (put_rels_in_dict(rels_dict,
                             rev_rels,
                             rev_read_pos,
                             rev_end3,
                             max(rev_end5, both_end3 + 1)))
            {return -1;}
    }

    // Fill relationships in the region of overlap.
    assert(both_end5 >= fwd_read_pos);
    assert(both_end5 >= rev_read_pos);
    unsigned char rel;
    for (size_t pos = both_end5; pos <= both_end3; pos++)
    {
        rel = fwd_rels[pos - fwd_read_pos] & rev_rels[pos - rev_read_pos];
        if (put_rel_in_dict(rels_dict, pos, rel))
            {return -1;}
    }

    return 0;
}


static int put_end_in_list(PyObject *ends_list, Py_ssize_t index, size_t end)
{
    PyObject *py_end = PyLong_FromSize_t(end);
    if (py_end == NULL) {return -1;}
    if (PyList_SetItem(ends_list, index, py_end)) {return -1;}
    return 0;
}


/*
Clean up by freeing all the dynamically allocated memory and optionally
decrementing the reference count for py_object, to prevent memory leaks.
*/
static PyObject *cleanup(unsigned char **rels1_ptr,
                         unsigned char **rels2_ptr,
                         IndelPodArray *pods,
                         PyObject *py_object)
{
    // All of the pointers to pointers must be non-NULL.
    assert(rels1_ptr != NULL);
    assert(rels2_ptr != NULL);
    assert(pods != NULL);

    // Free rels1 and point it to NULL.
    unsigned char *rels1 = *rels1_ptr;
    if (rels1 != NULL)
    {
        // printf("free %p\n", rels1);
        free(rels1);
        // Now that its memory has been freed, point rels1 at NULL.
        *rels1_ptr = NULL;
    }

    // Free rels2 and point it to NULL.
    unsigned char *rels2 = *rels2_ptr;
    if (rels2 != NULL)
    {
        // printf("free %p\n", rels2);
        free(rels2);
        // Now that its memory has been freed, point rels2 at NULL.
        *rels2_ptr = NULL;
    }

    // Free pods.
    free_pods(pods);

    // Decrement the reference count for py_object.
    if (py_object != NULL)
    {
        Py_DECREF(py_object);
    }

    // Return NULL so that py_calc_rels_lines can "return cleanup(...);"
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
    Py_ssize_t py_ref_len;
    unsigned long min_mapq;
    unsigned char min_qual;
    int insert3;
    int ambindel;
    int overhangs;
    unsigned long clip_end5;
    unsigned long clip_end3;
    if (!PyArg_ParseTuple(args,
                          "ssss#kbpppkk",
                          &line1,
                          &line2,
                          &ref,
                          &ref_seq,
                          &py_ref_len,
                          &min_mapq,
                          &min_qual,
                          &insert3,
                          &ambindel,
                          &overhangs,
                          &clip_end5,
                          &clip_end3))
        {return NULL;}
    size_t ref_len = (size_t)py_ref_len;
    
    // It is impossible for any of these pointers to be NULL.
    assert(line1 != NULL);
    assert(line2 != NULL);
    assert(ref != NULL);
    assert(ref_seq != NULL);

    // Determine the number of mates.
    Py_ssize_t num_mates = (*line2) ? 2 : 1;

    // Initialize containers to hold the results.
    
    unsigned char *rels1 = NULL, *rels2 = NULL;
    
    IndelPodArray pods;
    pods.capacity = 0;
    pods.pods = NULL;
    pods.num_pods = 0;
    
    PyObject *ends_rels_tuple = PyTuple_New(2);
    if (ends_rels_tuple == NULL)
        {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    
    PyObject *ends_tuple = PyTuple_New(2);
    if (ends_tuple == NULL)
        {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_rels_tuple, 0, ends_tuple);

    PyObject *end5s_list = PyList_New(num_mates);
    if (end5s_list == NULL)
        {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_tuple, 0, end5s_list);

    PyObject *end3s_list = PyList_New(num_mates);
    if (end3s_list == NULL)
        {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_tuple, 1, end3s_list);
    
    PyObject *rels_dict = PyDict_New();
    if (rels_dict == NULL)
        {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_rels_tuple, 1, rels_dict);

    // Calculate relationships for line 1.
    SamRead read1;
    if (calc_rels_line(&rels1,
                       &pods,
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
        {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}

    if (num_mates > 1)
    {
        // The read comprises forward (fwd) and reverse (rev) mates.
        size_t fwd_end5, fwd_end3, rev_end5, rev_end3;
        // Check if line 2 differs from line 1.
        if (strcmp(line1, line2))
        {
            // Free the pods for line 1 so that line 2 can reuse them.
            free_pods(&pods);

            // Calculate relationships for line 2.
            SamRead read2;
            if (calc_rels_line(&rels2,
                               &pods,
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
                {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}

            // Check if reads 1 and 2 are paired properly.
            if (validate_pair(&read1, &read2))
                {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}

            // Determine the 5'/3' ends of the forward/reverse reads.
            size_t fwd_read_pos, rev_read_pos;
            unsigned char *fwd_rels = NULL, *rev_rels = NULL;
            if (read2.reverse)
            {
                fwd_read_pos = read1.pos;
                fwd_end5 = read1.ref_end5;
                fwd_end3 = read1.ref_end3;
                fwd_rels = rels1;
                rev_read_pos = read2.pos;
                rev_end5 = read2.ref_end5;
                rev_end3 = read2.ref_end3;
                rev_rels = rels2;
            }
            else
            {
                fwd_read_pos = read2.pos;
                fwd_end5 = read2.ref_end5;
                fwd_end3 = read2.ref_end3;
                fwd_rels = rels2;
                rev_read_pos = read1.pos;
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
                                   fwd_read_pos,
                                   fwd_end5,
                                   fwd_end3,
                                   fwd_rels,
                                   rev_read_pos,
                                   rev_end5,
                                   rev_end3,
                                   rev_rels))
                {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
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
                                 read1.pos,
                                 read1.ref_end5,
                                 read1.ref_end3))
                {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
        }

        // Fill in the lists of 5' and 3' ends.
        if (put_end_in_list(end5s_list, 0, fwd_end5))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
        if (put_end_in_list(end5s_list, 1, rev_end5))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 0, fwd_end3))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 1, rev_end3))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    }
    else
    {
        // Line 2 does not exist.
        if (put_rels_in_dict(rels_dict,
                             rels1,
                             read1.pos,
                             read1.ref_end5,
                             read1.ref_end3))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
        
        // Fill int the lists of 5' and 3' ends.
        if (put_end_in_list(end5s_list, 0, read1.ref_end5))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 0, read1.ref_end3))
            {return cleanup(&rels1, &rels2, &pods, ends_rels_tuple);}
    }

    // Free the memory of rels1 and rels2, but not of the object that
    // will be returned.
    cleanup(&rels1, &rels2, &pods, NULL);

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
