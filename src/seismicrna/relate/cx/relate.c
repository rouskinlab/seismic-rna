// Required for Python <=3.12, before #include <Python.h>
#define PY_SSIZE_T_CLEAN
// Required for integration with Python; also includes other required
// header files such as <stdlib.h> and <string.h>.
#include <Python.h>




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Fundamental Definitions                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




static PyObject *RelateError;

static const int DECIMAL = 10;

static const char BASEA = 'A';
static const char BASEC = 'C';
static const char BASEG = 'G';
static const char BASET = 'T';


/* Check whether a character is a standard DNA base (A, C, G, or T). */
static inline int is_acgt(char base)
{
    return (base == BASEA || base == BASEC || base == BASEG || base == BASET);
}


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




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Relationship Encodings                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




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


/* Encode a base character as a substitution. */
static inline unsigned char encode_subs(char base)
{
    switch (base)
    {
        case BASEA:
            return SUB_A;
        case BASEC:
            return SUB_C;
        case BASEG:
            return SUB_G;
        case BASET:
            return SUB_T;
        default:
            // This function should only be called with A, C, G, or T.
            assert(0);
            return 0;
    }
}


/* Encode a match or substitution as such if the quality is sufficient,
otherwise as ambiguous. */
static inline unsigned char encode_relate(char ref_base,
                                          char read_base,
                                          char read_qual,
                                          unsigned char min_qual)
{
    // If the reference base is unknown, then the read base could be a
    // match or a substitution to any base.
    if (!is_acgt(ref_base))
        {return ANY_N;}
    // If the read base is unknown or low quality, then it could be a
    // match or substitution to any base but the reference base.
    if (read_qual < min_qual || !is_acgt(read_base))
        {return ANY_N ^ encode_subs(ref_base);}
    // Both reference and read bases are known and high-quality.
    return (read_base == ref_base) ? MATCH : encode_subs(read_base);
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Insertions and Deletions                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




static const int DELETION = 0;
static const int INSERTION = 1;
static const size_t INIT_POD_CAPACITY = 8;
static const size_t INIT_PODS_CAPACITY = 4;
static const size_t CAPACITY_FACTOR = 2;


static inline unsigned char get_ins_rel(int insert3)
{
    return insert3 ? INS_3 : INS_5;
}


static inline size_t calc_lateral5(size_t lateral3)
{
    assert(lateral3 >= 1);
    return lateral3 - 1;
}


static inline size_t get_lateral(size_t lateral3, int insert3)
{
    return insert3 ? lateral3 : calc_lateral5(lateral3);
}


typedef struct
{
    size_t label;
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


/* Initialize an IndelPod.
Do NOT call this function more than once on an IndelPod without first
calling free_pod(), or else memory in pod->indels will leak.
Make SURE to call free_pod() eventually, or else memory in pod->indels
will leak. */
static int init_pod(IndelPod *pod, size_t order, int insert)
{
    assert(pod != NULL);
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


/* Initialize an IndelPod.
Do NOT call this function on an IndelPod that has not been initialized
with init_pods(), or else the behavior is undefined and could free
random memory, possibly triggering an access violation or corrupting the
memory. */
static void free_pod(IndelPod *pod)
{
    assert(pod != NULL);
    free(pod->indels);
    // printf("free %p\n", pod->indels);
    pod->indels = NULL;
}


typedef struct
{
    IndelPod *pods;
    size_t capacity;
    size_t num_pods;
} IndelPodArray;


/* Initialize an IndelPodArray.
Do NOT call this function on an IndelPodArray after calling add_indel(),
unless free_pods() is called after add_indel() and before init_pods(),
or else memory in pods->pods will leak. */
static void init_pods(IndelPodArray *pods)
{
    assert(pods != NULL);
    // Initialize pods->pods to NULL so that the pointer doesn't dangle.
    pods->pods = NULL;
    pods->capacity = 0;
    pods->num_pods = 0;
}


/* Free the dynamically allocated memory for an IndelPodArray and all
IndelPods in the array.
Do NOT call this function on an IndelPodArray that has not yet been
initialized with init_pods(), or else the behavior is undefined and
could free random memory, possibly triggering an access violation or
corrupting the memory. */
static void free_pods(IndelPodArray *pods)
{
    assert(pods != NULL);
    // Free each IndelPod.
    for (size_t p = 0; p < pods->num_pods; p++)
    {
        free_pod(&(pods->pods[p]));
    }
    free(pods->pods);
    // printf("free %p\n", pods->pods);
    // Now that its memory has been freed, reset all of its fields.
    init_pods(pods);
}


/* Add an Indel to the last IndelPod in an IndelPodArray. If the array
has no IndelPods, or if the last IndelPod has the wrong type of Indel,
then automatically append a new IndelPod to the IndelPodArray.
Do NOT call this function on an IndelPodArray that has not yet been
initialized with init_pods(), or else the behavior is undefined and
could free random memory, possibly triggering an access violation or
corrupting the memory.
Make SURE to call free_pods() eventually, or else memory in pods->pods
will leak. */
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
            // The last pod is the wrong type: a new pod is needed.
            if (pods->num_pods == pods->capacity)
            {
                // The array of pods is already at its maximum capacity:
                // allocate a new array with twice the capacity.
                pods->capacity *= CAPACITY_FACTOR;
                IndelPod *new_pods = realloc(pods->pods,
                                             pods->capacity * sizeof(IndelPod));
                if (new_pods == NULL)
                {
                    PyErr_NoMemory();
                    return -1;
                }
                // printf("implicit-free %p\n", pods->pods);
                // printf("realloc %p\n", new_pods);
                pods->pods = new_pods;
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
                Indel *new_indels = realloc(pod->indels,
                                            pod->capacity * sizeof(Indel));
                if (new_indels == NULL)
                {
                    PyErr_NoMemory();
                    return -1;
                }
                // printf("implicit-free %p\n", pod->indels);
                // printf("realloc %p\n", new_indels);
                pod->indels = new_indels;
            }
        }
    }

    // Initialize the new indel.
    assert(pod != NULL);
    assert(pod->num_indels < pod->capacity);
    indel = &(pod->indels[pod->num_indels]);
    indel->label = pod->num_indels;
    indel->insert = insert;
    indel->opposite = opposite;
    indel->lateral3 = lateral3;
    pod->num_indels++;
    return 0;
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// SAM Line Parser                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




static const char *SAM_SEP = "\t\n";
static const char CIG_ALIGN = 'M';
static const char CIG_DELET = 'D';
static const char CIG_INSRT = 'I';
static const char CIG_MATCH = '=';
static const char CIG_SCLIP = 'S';
static const char CIG_SUBST = 'X';
static const unsigned long FLAG_PAIRED = 1;
static const unsigned long FLAG_PROPER = 2;
static const unsigned long FLAG_REV = 16;
static const unsigned long FLAG_READ1 = 64;
static const unsigned long FLAG_READ2 = 128;
static const unsigned long MAX_FLAG = 4095;  // 4095 = 2^12 - 1


typedef struct
{
    const char *op;  // type of CIGAR operation
    size_t len;      // length of CIGAR operation
} CigarOp;


/* Initialize a CigarOp. */
static inline void init_cigar(CigarOp *cigar, const char *cigar_str)
{
    assert(cigar_str != NULL);
    // Point the operation to one position before the CIGAR string.
    cigar->op = cigar_str - 1;
    cigar->len = 0;
}


/* Find the next operation in a CIGAR string and store it in cigar.
Do NOT call this function on a CigarOp that has not yet been initialized
with init_cigar(), or else it will try to read random memory, generating
incorrect results and possibly causing an access violation. */
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
            PyErr_SetString(RelateError, "Invalid CIGAR operation");
            return -1;
        }
        // Point cigar->op to NULL to signify reaching the end.
        cigar->op = NULL;
    }
    return 0;
}


typedef struct
{
    // Source attributes
    char *line;           // buffer to hold line from SAM file
    // Parsed attributes
    char *name;           // read name
    unsigned long flag;   // bitwise flag
    char *ref;            // reference name
    size_t pos;           // mapping position in the ref (1-indexed)
    unsigned long mapq;   // mapping quality
    char *cigar;          // CIGAR string
    char *seq;            // read sequence
    char *quals;          // read quality scores
    // Calculated attributes
    char *end;            // pointer to address after 3' end of read
    size_t len;           // read length
    int paired;           // whether read is paired
    int proper;           // whether read is properly paired
    int reverse;          // whether read is reverse complemented
    int read1;            // whether read is the 1st read
    int read2;            // whether read is the 2nd read
    size_t num_rels;      // number of relationships
    size_t ref_end5;      // 5' end of the read in the ref (1-indexed)
    size_t ref_end3;      // 3' end of the read in the ref (1-indexed)
    // Container attributes
    unsigned char *rels;  // relationship array
    IndelPodArray pods;   // pods for the ambiguous indel algorithm
} SamRead;


/* Initialize the fields of a SamRead using a line from a SAM file.
Do NOT call this function more than once on the same SamRead without
first calling free_read(), or else the memory allocated for read->rels
and read->pods will leak. */
static void init_read(SamRead *read)
{
    assert(read != NULL);
    // Initialize all pointer fields to NULL to prevent errors caused
    // by uninitialized pointers.
    read->line = NULL;
    read->name = NULL;
    read->ref = NULL;
    read->cigar = NULL;
    read->seq = NULL;
    read->quals = NULL;
    read->end = NULL;
    read->rels = NULL;
    // Initialize the read's pods.
    init_pods(&(read->pods));
}


/* Free all of the dynamically allocated memory for a SamRead.
Do NOT call this function on a SamRead that has not yet been initialized
with init_read(), or else read->line and read->rels will be freed before
being initialized, possibly causing an access violation or corrupting
the memory. */
static void free_read(SamRead *read)
{
    assert(read != NULL);
    free(read->line);
    // printf("free %p\n", read->line);
    read->line = NULL;
    free(read->rels);
    // printf("free %p\n", read->rels);
    read->rels = NULL;
    free_pods(&(read->pods));
    // None of the other fields in the SamRead are dynamically allocated
    // memory, so nothing else needs to be done with them.
}


/* Return the next SAM field in the string, or NULL on failure. */
static char *next_sam_field(char **end)
{
    assert(end != NULL);
    if (*end == NULL)
    {
        // The previous field was the last one possible, hence this one
        // is invalid.
        return NULL;
    }
    if (**end == *SAM_SEP)
    {
        // The previous field ended on the field separator, which means
        // that the current field is an empty string, which is invalid.
        return NULL;
    }
    // Return the beginning of the next field and advance *end to the 
    // end of the next field.
    // If the return value is NULL, it indicates that strtok_r() failed
    // to find the start and end positions of the next field, which is
    // an error (but this function already signals failure with NULL).
    // Note that *end is now allowed to be NULL, which means that this
    // field is the last in the string (which is valid in this case);
    // however, calling next_sam_field() on the same string one more
    // time will produce an error because there will be no next field.
    
    return strtok_r(NULL, SAM_SEP, end);
}


/* Parse str (null-terminated) into an unsigned integer ulong. */
static int parse_ulong(unsigned long *ulong, const char *str)
{
    assert(ulong != NULL);
    if (str == NULL)
    {
        // This conditional is included so that this function can be
        // called using parse_ulong(&ulong, parse_sam_field(&end));
        // parse_sam_field returns NULL if it fails, so this function
        // needs to handle (str == NULL) by propagating that error.
        PyErr_SetString(RelateError, "(Got NULL string to parse as ulong)");
        return -1;
    }

    // Parse str as a base-10 integer; store its numeric value in ulong.
    char *endptr;
    *ulong = strtoul(str, &endptr, DECIMAL);
    assert(endptr != NULL);

    // Check if the parse processed any characters.
    if (endptr == str)
    {
        assert(*ulong == 0);  // strtoul returns 0 on failure
        PyErr_SetString(RelateError, "(Parsed zero characters as ulong)");
        return -1;
    }
    // Check if all of the field was parsed; if so, then endptr will
    // point to the null terminator ('\0') at the end of str.
    if (*endptr)
    {
        PyErr_SetString(RelateError,
                        "(Failed to parse all characters as ulong)");
        return -1;
    }

    return 0;
}


/* Fill the fields of a SamRead using a line from a SAM file.
Do NOT call this function on a SamRead that has not yet been initialized
with init_read(), so that in case this function fails and returns -1,
the SamRead's memory can be freed safely with free_read(). */
static int parse_sam_line(SamRead *read,
                          const char *line,
                          size_t line_length,
                          size_t clip_end5,
                          size_t clip_end3,
                          size_t ref_len)
{
    assert(line != NULL);
    assert(read != NULL);
    // If any of these are != NULL, then memory can leak.
    assert(read->line == NULL);
    assert(read->rels == NULL);
    assert(read->pods.pods == NULL);
    // Passing line_length as an argument avoids needing to calculate it
    // with strlen, but will cause bugs if line_length != strlen(line).
    assert(strlen(line) == line_length);

    // Allocate (line_length + 1) chars to include the null terminator.
    size_t line_size = (line_length + 1) * sizeof(char);
    read->line = malloc(line_size);
    if (read->line == NULL)
    {
        PyErr_NoMemory();
        return -1;
    }
    // printf("malloc %p\n", read->line);

    // Copy the line into read->line and ensure it was copied correctly.
    memcpy(read->line, line, line_size);
    assert(read->line[line_length] == '\0');
    assert(strlen(read->line) == line_length);

    char *end = read->line;  // End of the current field.
    unsigned long temp_ulong;  // Hold parsed numbers before casting.

    // Read name
    if ((read->name = next_sam_field(&end)) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse read name");
        return -1;
    }

    // Bitwise flag
    if (parse_ulong(&temp_ulong, next_sam_field(&end)))
    {
        PyErr_SetString(RelateError, "Failed to parse SAM flag");
        return -1;
    }
    read->flag = temp_ulong;
    if (read->flag > MAX_FLAG)
    {
        PyErr_SetString(RelateError, "SAM flag is too large");
        return -1;
    }
    // Individual flag bits
    read->paired = (read->flag & FLAG_PAIRED) != 0;
    read->proper = (read->flag & FLAG_PROPER) != 0;
    read->reverse = (read->flag & FLAG_REV) != 0;
    read->read1 = (read->flag & FLAG_READ1) != 0;
    read->read2 = (read->flag & FLAG_READ2) != 0;

    // Reference name
    if ((read->ref = next_sam_field(&end)) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse reference name");
        return -1;
    }

    // Mapping position
    if (parse_ulong(&temp_ulong, next_sam_field(&end)))
    {
        PyErr_SetString(RelateError, "Failed to parse mapping position");
        return -1;
    }
    read->pos = (size_t)temp_ulong;
    if (read->pos == 0)
    {
        PyErr_SetString(RelateError, "Mapping position is 0");
        return -1;
    }
    if (read->pos > ref_len)
    {
        PyErr_SetString(RelateError,
                        "Mapping position exceeds length of reference");
        return -1;
    }

    // Mapping quality
    if (parse_ulong(&temp_ulong, next_sam_field(&end)))
    {
        PyErr_SetString(RelateError, "Failed to parse mapping quality");
        return -1;
    }
    read->mapq = temp_ulong;

    // CIGAR string
    if ((read->cigar = next_sam_field(&end)) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse CIGAR string");
        return -1;
    }

    // Next reference (ignored)
    if (next_sam_field(&end) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse next reference name");
        return -1;
    }

    // Next position (ignored)
    if (next_sam_field(&end) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse next mapping position");
        return -1;
    }

    // Template length (ignored)
    if (next_sam_field(&end) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse template length");
        return -1;
    }

    // Read sequence
    if ((read->seq = next_sam_field(&end)) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse read sequence");
        return -1;
    }

    // Read end position
    // Subtract 1 because end points to one after the \t delimiter, but
    // read->end must point to the tab delimiter itself, i.e.
    // ACGT FFFF      (length = 4)
    // ^read->seq
    //     ^read->end (= 4 + read->seq)
    //      ^end      (= 5 + read->seq)
    // It is possible here for end to be NULL, if the line ended at the
    // read sequence (with the quality string missing). If so, (end - 1)
    // will overflow, but that's okay because this function will catch
    // the error and return -1 when it tries to parse the read quality.
    read->end = end - 1;

    // Read length (using pointer arithmetic)
    read->len = read->end - read->seq;

    // Read quality
    if ((read->quals = next_sam_field(&end)) == NULL)
    {
        PyErr_SetString(RelateError, "Failed to parse read quality");
        return -1;
    }

    // Lengths of read and quality strings must match.
    // Subtract 1 from end for the same reason as in "Read end position"
    if (end != NULL && read->len != (size_t)((end - 1) - read->quals))
    {
        PyErr_SetString(RelateError,
                        "Read sequence and quality strings differ in length");
        return -1;
    }

    // Number of reference and read bases consumed by the CIGAR string.
    read->num_rels = 0;
    size_t num_read_bases = 0;
    // Initialize the CIGAR operation.
    CigarOp cigar;
    init_cigar(&cigar, read->cigar);
    // Read the first operation from the CIGAR string; catch errors.
    if (get_next_cigar_op(&cigar))
        {return -1;}
    // Read the entire CIGAR string one operation at a time.
    while (cigar.op != NULL)
    {
        char op = *cigar.op;
        // Decide what to do based on the current CIGAR operation.
        switch (op)
        {
            case CIG_ALIGN:
            case CIG_MATCH:
            case CIG_SUBST:
                // These operations consume the reference and read.
                read->num_rels += cigar.len;
                num_read_bases += cigar.len;
                break;
            case CIG_DELET:
                // These operations consume only the reference.
                read->num_rels += cigar.len;
                break;
            case CIG_INSRT:
            case CIG_SCLIP:
                // These operations consume only the read.
                num_read_bases += cigar.len;
                break;
            default:
                // Other operations are not supported.
                PyErr_SetString(RelateError,
                                "Unsupported CIGAR operation");
                return -1;
        }
        // Read the next CIGAR operation; catch errors.
        if (get_next_cigar_op(&cigar))
            {return -1;}
        if (cigar.op != NULL)
        {
            if (*cigar.op == op)
            {
                PyErr_SetString(RelateError,
                                "Identical consecutive CIGAR operations");
                return -1;
            }
            if ((*cigar.op == CIG_DELET && op == CIG_INSRT)
                ||
                (*cigar.op == CIG_INSRT && op == CIG_DELET))
            {
                PyErr_SetString(RelateError,
                                "Adjacent insertion and deletion");
                return -1;
            }
        }
    }
    // The number of reference bases consumed must be ≥ 1 but not extend
    // out of the reference sequence.
    if (read->num_rels == 0)
    {
        PyErr_SetString(RelateError,
                        "CIGAR operations consumed 0 bases in the reference");
        return -1;
    }
    if (read->num_rels + (read->pos - 1) > ref_len)
    {
        PyErr_SetString(RelateError,
                        "CIGAR operations extended out of the reference");
        return -1;
    }
    // The number of read bases consumed must equal the read length.
    if (num_read_bases != read->len)
    {
        PyErr_SetString(RelateError,
                        "CIGAR operations consumed a number of read bases "
                        "different from the read length");
        return -1;
    }

    // Reference position at which the 5' end of the read is located,
    // which is the mapping position plus the 5' clip (up to a maximum
    // of the reference length plus 1).
    read->ref_end5 = min(read->pos + clip_end5, ref_len + 1);

    // Reference position at which the 3' end of the read is located,
    // which is the mapping position plus one less than the number of
    // reference bases consumed, minus the 3' clip (down to a minimum
    // of 0). Because read->num_rels must be ≥ 1, subtracting 1 will
    // not yield a negative number and cause overflow.
    assert(read->num_rels >= 1);
    size_t ref_end3 = read->pos + (read->num_rels - 1);
    // Also to avoid overflow, subtract clip_end3 only after confirming
    // that the difference would not be negative.
    read->ref_end3 = (ref_end3 > clip_end3) ? (ref_end3 - clip_end3) : 0;

    return 0;
}


/* Validate that a SamRead's reference name matches the name of the SAM
file and that its mapping quality is sufficient. */
static int validate_read(const SamRead *read,
                         const char *ref,
                         unsigned long min_mapq,
                         int paired,
                         int proper)
{
    assert(read != NULL);
    assert(read->ref != NULL);
    assert(ref != NULL);

    // Validate alignment attributes.
    if (strcmp(read->ref, ref))
    {
        PyErr_SetString(RelateError,
                        "Reference name does not match name of SAM file");
        return -1;
    }
    if (read->mapq < min_mapq)
    {
        PyErr_SetString(RelateError,
                        "Mapping quality is insufficient");
        return -1;
    }

    // Validate paired-end attributes.
    if (read->paired != paired)
    {
        if (paired)
        {
            PyErr_SetString(RelateError,
                            "Lines indicate read should be paired-end, "
                            "but it is marked as single-end");
            return -1;
        }
        else
        {
            PyErr_SetString(RelateError,
                            "Lines indicate read should be single-end, "
                            "but it is marked as paired-end");
            return -1;
        }
    }
    if (read->proper != proper)
    {
        if (proper)
        {
            PyErr_SetString(RelateError,
                            "Lines indicate read should be properly paired, "
                            "but it is marked as improperly paired");
            return -1;
        }
        else
        {
            PyErr_SetString(RelateError,
                            "Lines indicate read should be improperly paired, "
                            "but it is marked as properly paired");
            return -1;
        }
    }

    return 0;
}


/* Validate that two paired-end SamReads that are allegedly mates have
the same name and other compatible attributes. */
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
        PyErr_SetString(RelateError, "Mates 1 and 2 have different names");
        return -1;
    }

    // Mates 1 and 2 must be marked as READ1 and READ2, respectively.
    if (!(read1->read1) || read1->read2)
    {
        PyErr_SetString(RelateError, "Mate 1 is not marked as READ1");
        return -1;
    }
    if (!(read2->read2) || read2->read1)
    {
        PyErr_SetString(RelateError, "Mate 2 is not marked as READ2");
        return -1;
    }

    // Mates 1 and 2 must be in opposite orientations (i.e. one must be
    // marked as REVERSE and the other as not REVERSE).
    if (read1->reverse == read2->reverse)
    {
        PyErr_SetString(RelateError,
                        "Mates 1 and 2 aligned in the same orientation");
        return -1;
    }

    return 0;
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Ambiguous Insertions and Deletions Algorithm                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




/* Determine the order of two IndelPods from 5' to 3'. */
static inline int comp_pods_5to3(const void *a, const void *b)
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


/* Determine the order of two IndelPods from 3' to 5'. */
static inline int comp_pods_3to5(const void *a, const void *b)
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
static inline void sort_pods(IndelPodArray *pods, int move5to3)
{
    assert(pods != NULL);
    if (move5to3)
        {qsort(pods->pods, pods->num_pods, sizeof(IndelPod), comp_pods_5to3);}
    else
        {qsort(pods->pods, pods->num_pods, sizeof(IndelPod), comp_pods_3to5);}
}


/* Determine the order of two Indels from 5' to 3'. */
static inline int comp_indels_5to3(const void *a, const void *b)
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


/* Determine the order of two Indels from 3' to 5'. */
static inline int comp_indels_3to5(const void *a, const void *b)
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
static inline void sort_pod(IndelPod *pod, int move5to3)
{
    assert(pod != NULL);
    if (move5to3)
        {qsort(pod->indels, pod->num_indels, sizeof(Indel), comp_indels_5to3);}
    else
        {qsort(pod->indels, pod->num_indels, sizeof(Indel), comp_indels_3to5);}
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
                        size_t label,
                        size_t opposite,
                        size_t lateral3)
{
    // Locate the indel with the given label.
    Indel *indel = NULL;
    size_t i = 0;
    while (indel == NULL && i < pod->num_indels)
    {
        if (pod->indels[i].label == label)
        {
            indel = &(pod->indels[i]);
        }
        i++;
    }
    assert(indel != NULL);
    assert(indel->insert == pod->insert);
    assert(indel->opposite != opposite);

    // Move every other indel in the pod that lies between the given
    // indel and the position to which the given indel should be moved.
    for (i = 0; i < pod->num_indels; i++)
    {
        Indel *other = &(pod->indels[i]);
        if (other != indel)
        {
            assert(other->insert == pod->insert);
            assert(other->label != label);
            // Check if the other indel lies between the given indel and
            // the position to which the given indel should be moved.
            assert(other->opposite != indel->opposite);
            assert(other->opposite != opposite);
            if ((indel->opposite < other->opposite)
                ==
                (other->opposite < opposite))
            {
                // If so, then move it.
                move_indel(other, other->opposite, lateral3);
            }
        }
    }

    // Move the given indel.
    move_indel(indel, opposite, lateral3);
}


/* Return whether a pod has an indel with the given opposite value. */
static inline int pod_has_opposite(const IndelPod *pod, size_t opposite)
{
    for (size_t i = 0; i < pod->num_indels; i++)
    {
        if (pod->indels[i].opposite == opposite)
            {return 1;}
    }
    return 0;
}


static void calc_positions(size_t *swap_lateral,
                           size_t *next_lateral3,
                           size_t *next_opposite,
                           const IndelPod *pod,
                           const Indel *indel,
                           int move5to3)
{
    assert(swap_lateral != NULL);
    assert(next_lateral3 != NULL);
    assert(next_opposite != NULL);
    assert(pod != NULL);
    assert(pod->indels != NULL);
    assert(indel != NULL);
    assert(indel >= pod->indels
           && (size_t)(indel - pod->indels) < pod->num_indels);

    // Find the position within the indel's own sequence with which to
    // try to swap the indel.
    *swap_lateral = get_lateral(indel->lateral3, move5to3);

    // Determine the direction in which to move the indel.
    int direction = move5to3 ? 1 : -1;

    if (direction == -1)
    {
        // These attributes must be ≥ 1 or else they will overflow.
        assert(indel->opposite >= 1);
        assert(indel->lateral3 >= 1);
    }

    // If the indel moves, then its position within its own sequence
    // will change.
    *next_lateral3 = indel->lateral3 + direction;

    // Find the position within the opposite sequence to which to try
    // to move the indel.
    *next_opposite = indel->opposite + direction;
    while (pod_has_opposite(pod, *next_opposite))
    {
        if (direction == -1)
        {
            // Must be ≥ 1 or else it will overflow.
            assert(*next_opposite >= 1);
        }
        *next_opposite += direction;
    }
}


/* Check whether the indel would collide with another pod. */
static inline int check_collisions(const IndelPodArray *pods,
                                   size_t pod_index,
                                   size_t next_lateral3)
{
    size_t next_pod_index = pod_index + 1;
    if (next_pod_index < pods->num_pods)
    {
        IndelPod *next_pod = &(pods->pods[next_pod_index]);
        // Confirm the next pod has a different type of indel.
        assert(next_pod->insert != pods->pods[pod_index].insert);
        // Every pod must contain at least 1 indel.
        assert(next_pod->num_indels >= 1);
        // In the next pod, check the first indel, which is the one that
        // is closest to the current indel and could collide with it.
        Indel *next_indel = &(next_pod->indels[0]);
        return (next_indel->opposite == next_lateral3
                ||
                next_indel->opposite == calc_lateral5(next_lateral3));
    }
    // There is no next pod with which to collide.
    return 0;
}


/* Check if two relationships are consistent with each other. */
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


/* Calculate relationships for a deletion. */
static void calc_rels_del(unsigned char *rel_del,
                          unsigned char *rel_opp,
                          const char *rels_seq,
                          const char *read_seq,
                          const char *read_quals,
                          unsigned char min_qual,
                          size_t swap_lateral,
                          size_t next_opposite,
                          size_t opposite)
{
    assert(rel_del != NULL);
    assert(rel_opp != NULL);
    assert(rels_seq != NULL);
    assert(read_seq != NULL);
    assert(read_quals != NULL);
    // Read base and its quality with which the deletion will swap.
    char read_base = read_seq[swap_lateral];
    char read_qual = read_quals[swap_lateral];
    // Relationship for the base that is currently deleted from the 
    // read, which would move opposite the read base.
    *rel_del = encode_relate(rels_seq[opposite],
                             read_base,
                             read_qual,
                             min_qual);
    // Relationship for the base that is currently opposite the read
    // base, which would become a deletion after the move.
    *rel_opp = encode_relate(rels_seq[next_opposite],
                             read_base,
                             read_qual,
                             min_qual);
}


/* Calculate relationships for a insertion. */
static void calc_rels_ins(unsigned char *rel_ins,
                          unsigned char *rel_opp,
                          const char *rels_seq,
                          const char *read_seq,
                          const char *read_quals,
                          unsigned char min_qual,
                          size_t swap_lateral,
                          size_t next_opposite,
                          size_t opposite)
{
    assert(rel_ins != NULL);
    assert(rel_opp != NULL);
    assert(rels_seq != NULL);
    assert(read_seq != NULL);
    assert(read_quals != NULL);
    // Reference base with which the insertion will swap.
    char ref_base = rels_seq[swap_lateral];
    // Relationship for the base that is currently inserted into the
    // read, which would move opposite the reference base.
    *rel_ins = encode_relate(ref_base,
                             read_seq[opposite],
                             read_quals[opposite],
                             min_qual);
    // Relationship for the base that is currently opposite the 
    // reference base, which would become an insertion after the move.
    *rel_opp = encode_relate(ref_base,
                             read_seq[next_opposite],
                             read_quals[next_opposite],
                             min_qual);
}


/* Determine if the position adjacent to an insertion must be marked. */
static int need_mark_ins_adj(unsigned char rel,
                             const IndelPod *pod,
                             size_t init_lateral,
                             int insert3)
{
    assert(pod != NULL);
    // The position must be marked if the relationship is not a match.
    if (rel != MATCH)
        {return 1;}
    // If it is a match, then it must be marked as such unless any
    // insertion is still adjacent to it on the side opposite the
    // insertion that just moved (because then it would count as an
    // insertion, not as a match).
    for (size_t i = 0; i < pod->num_indels; i++)
    {
        if (get_lateral(pod->indels[i].lateral3, insert3) == init_lateral)
            {return 0;}
    }
    return 1;
}


/* Try to move an indel. */
static int try_move_indel(Indel *indel,
                          SamRead *read,
                          const char *rels_seq,
                          size_t read_end5,
                          size_t read_end3,
                          unsigned char min_qual,
                          int insert3,
                          int move5to3,
                          size_t pod_index)
{
    // Select the pod.
    assert(read != NULL);
    assert(pod_index < read->pods.num_pods);
    IndelPod *pod = &(read->pods.pods[pod_index]);
    assert(indel != NULL
           && indel >= pod->indels
           && (size_t)(indel - pod->indels) < pod->num_indels);
    assert(insert3 == 0 || insert3 == 1);
    assert(move5to3 == 0 || move5to3 == 1);
    
    // Calculate the positions to which to move the indel.
    size_t swap_lateral, next_lateral3, next_opposite;
    calc_positions(&swap_lateral,
                   &next_lateral3,
                   &next_opposite,
                   pod,
                   indel,
                   move5to3);

    // Choose the method based on the type of indel.
    if (indel->insert)
    {
        // The indel is an insertion.
        // Stop if the next position would be the first or last position
        // in the non-soft-clipped region of the read.
        // Both read_end5 and read_end3 are 1-indexed, but next_opposite
        // is 0-indexed.
        assert(read_end5 >= 1);
        assert(read_end3 >= 1);
        assert(read_end3 <= read->len);
        if (next_opposite < read_end5 || next_opposite >= read_end3 - 1)
            {return 0;}
        // Stop if the indel would collide with another pod.
        if (check_collisions(&(read->pods),
                             pod_index,
                             next_lateral3))
            {return 0;}
        
        // Calculate what the relationships would be after the move.
        unsigned char rel_indel, rel_opp;
        calc_rels_ins(&rel_indel,
                      &rel_opp,
                      rels_seq,
                      read->seq,
                      read->quals,
                      min_qual,
                      swap_lateral,
                      next_opposite,
                      indel->opposite);
        // Stop if the relationships before and after the move are not
        // consistent with each other.
        if (!consistent_rels(rel_indel, rel_opp))
            {return 0;}
        
        // Initial lateral position of the insertion.
        size_t init_lateral = get_lateral(indel->lateral3, insert3);

        // Relationship code for insertions.
        unsigned char ins_rel = get_ins_rel(insert3);

        if (move5to3 == insert3)
        {
            // The insertion moves to the same side as it is marked on.
            // Move the insertion.
            move_indels(pod, indel->label, next_opposite, next_lateral3);
            assert(get_lateral(indel->lateral3, insert3) != init_lateral);
            
            if (need_mark_ins_adj(rel_indel, pod, init_lateral, insert3))
            {
                // Mark the position with which the insertion swapped.
                read->rels[get_lateral(indel->lateral3, !insert3)] |= rel_indel;
            }

            // Mark the position to which the insertion moved.
            read->rels[get_lateral(indel->lateral3, insert3)] |= ins_rel;
        }
        else
        {
            // The insertion moves to the side opposite that on which it
            // is marked.
            // Determine the relationship to mark at the position that
            // is adjacent to the insertion on the side opposite the
            // direction of movement.
            size_t swap_lateral_, next_lateral3_, next_opposite_;
            calc_positions(&swap_lateral_,
                           &next_lateral3_,
                           &next_opposite_,
                           pod,
                           indel,
                           insert3);
            unsigned char rel_indel_, rel_opp_;
            calc_rels_ins(&rel_indel_,
                          &rel_opp_,
                          rels_seq,
                          read->seq,
                          read->quals,
                          min_qual,
                          swap_lateral_,
                          next_opposite_,
                          indel->opposite);
            
            // Move the insertion.
            move_indels(pod, indel->label, next_opposite, next_lateral3);
            assert(get_lateral(indel->lateral3, insert3) != init_lateral);

            if (need_mark_ins_adj(rel_opp_, pod, init_lateral, insert3))
            {
                // Mark the position that was adjacent to the insertion
                // before it moved.
                read->rels[init_lateral] |= rel_opp_;
            }

            // Mark the position to which the insertion moved.
            size_t curr_lateral = get_lateral(indel->lateral3, insert3);
            read->rels[curr_lateral] |= ins_rel;
            if (rel_indel != MATCH)
            {
                read->rels[curr_lateral] |= rel_indel;
            }
        }
    }
    else
    {
        // The indel is a deletion.
        // Stop if the next position would be the first or last position
        // in the reference.
        assert(read->num_rels >= 1);
        if (next_opposite == 0 || next_opposite >= read->num_rels - 1)
            {return 0;}
        // Stop if the indel would collide with another pod.
        if (check_collisions(&(read->pods),
                             pod_index,
                             next_lateral3))
            {return 0;}
        
        // Calculate what the relationships would be after the move.
        unsigned char rel_indel, rel_opp;
        calc_rels_del(&rel_indel,
                      &rel_opp,
                      rels_seq,
                      read->seq,
                      read->quals,
                      min_qual,
                      swap_lateral,
                      next_opposite,
                      indel->opposite);
        // Stop if the relationships before and after the move are not
        // consistent with each other.
        if (!consistent_rels(rel_indel, rel_opp))
            {return 0;}
        
        // When the deletion moves, the position from which it moves
        // (indel->opposite) gains the relationship (rel_del) between
        // the read base and the reference base that was originally
        // deleted from the read.
        read->rels[indel->opposite] |= rel_indel;
        // Move the deletion.
        move_indels(pod, indel->label, next_opposite, next_lateral3);
        // Mark the position to which the deletion moved.
        read->rels[next_opposite] |= DELET;
    }

    // Re-sort the indels in the pod so they stay in positional
    // order following the move.
    sort_pod(pod, move5to3);

    // Signify that the indel moved.
    return 1;
}


/* Recursive algorithm to mark all ambiguous indels. */
static void find_ambindels_recurse(SamRead *read,
                                   const char *rels_seq,
                                   size_t read_end5,
                                   size_t read_end3,
                                   unsigned char min_qual,
                                   int insert3,
                                   int move5to3,
                                   size_t pod_index,
                                   size_t indel_index)
{
    // Select the pod at this index.
    assert(pod_index < read->pods.num_pods);
    IndelPod *pod = &(read->pods.pods[pod_index]);
    // Select the indel at this index.
    assert(indel_index < pod->num_indels);
    Indel *indel = &(pod->indels[indel_index]);

    // Record the initial attributes of the indel.
    size_t label = indel->label;
    size_t init_opposite = indel->opposite;
    size_t init_lateral3 = indel->lateral3;

    // Try to move the indel.
    if (try_move_indel(indel,
                       read,
                       rels_seq,
                       read_end5,
                       read_end3,
                       min_qual,
                       insert3,
                       move5to3,
                       pod_index))
    {
        // Try to move the indel at this index another step.
        find_ambindels_recurse(read,
                               rels_seq,
                               read_end5,
                               read_end3,
                               min_qual,
                               insert3,
                               move5to3,
                               pod_index,
                               indel_index);
        
        // Backtracking is only needed when moving from 5' to 3'.
        if (move5to3)
        {
            // Move the indel back to its initial position.
            move_indels(pod, label, init_opposite, init_lateral3);
            // Re-sort the indels in the pod so they stay in positional
            // order after moving back to their initial positions.
            sort_pod(pod, move5to3);
        }
    }

    // Check if the pod contains another indel after the current one.
    if (indel_index + 1 < pod->num_indels)
    {
        // Move the next indel in the pod.
        find_ambindels_recurse(read,
                               rels_seq,
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
        find_ambindels_recurse(read,
                               rels_seq,
                               read_end5,
                               read_end3,
                               min_qual,
                               insert3,
                               move5to3,
                               pod_index - 1,
                               0);
    }
}


/* Find and mark all ambiguous insertions and deletions. */
static void find_ambindels(SamRead *read,
                           const char *rels_seq,
                           size_t read_end5,
                           size_t read_end3,
                           unsigned char min_qual,
                           int insert3)
{
    if (read->pods.num_pods == 0)
    {
        // Nothing to do.
        return;
    }
    assert(read->pods.pods != NULL);
    
    // This algorithm works in two stages:
    // 1. Move every indel as far as possible in the 5' direction
    //    (move5to3 == 0) and leave them there for the next stage.
    // 2. Move every indel as far as possible in the 3' direction
    //    (move5to3 == 1) using backtracking to make sure that all
    //    possible movements are considered.
    for (int move5to3 = 0; move5to3 <= 1; move5to3++)
    {
        // Sort all of the pods.
        sort_pods(&(read->pods), move5to3);
        // In each pod, sort the indels.
        for (size_t p = 0; p < read->pods.num_pods; p++)
        {
            sort_pod(&(read->pods.pods[p]), move5to3);
        }
        find_ambindels_recurse(read,
                               rels_seq,
                               read_end5,
                               read_end3,
                               min_qual,
                               insert3,
                               move5to3,
                               read->pods.num_pods - 1,
                               0);
    }
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Relationship Calculator                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




/* Compute the relationships for a line in a SAM file.
Do NOT call this function on a SamRead that has not yet been initialized
with init_read(), or else read->pods will not be initialized properly,
which will cause undefined and probably erroneous behavior.
Do NOT call this function more than once on the same SamRead without
first calling free_read(), or else the memory originally allocated to
read->rels and read->pods will leak.
Make SURE to call free_read() eventually, or else memory in pods->rels
will leak. */
static int calc_rels_line(SamRead *read,
                          const char *line,
                          size_t line_length,
                          const char *ref,
                          const char *ref_seq,
                          size_t ref_len,
                          unsigned long min_mapq,
                          int paired,
                          int proper,
                          unsigned char min_qual,
                          int insert3,
                          int ambindel,
                          size_t clip_end5,
                          size_t clip_end3)
{
    // Validate the SAM file line and parse it into a SamRead.
    assert(read != NULL);
    assert(line != NULL);
    if (parse_sam_line(read, line, line_length, clip_end5, clip_end3, ref_len))
        {return -1;}
    if (validate_read(read, ref, min_mapq, paired, proper))
        {return -1;}

    // Point to where the read starts in the reference sequence.
    // Because read->pos is guaranteed to be ≥ 1, subtracting 1 will not
    // yield a negative number and cause overflow.
    assert(read->pos >= 1);
    const char *rels_seq = ref_seq + (read->pos - 1);

    // Allocate memory for the relationships. The relationships will be
    // filled using the CIGAR string, so the number of relationships
    // allocated must equal the number of reference bases consumed by
    // the CIGAR string, not the number after clipping.
    assert(read->rels == NULL);  // If != NULL, then memory will leak.
    assert(read->num_rels >= 1);
    read->rels = malloc(read->num_rels * sizeof(unsigned char));
    if (read->rels == NULL)
    {
        PyErr_NoMemory();
        return -1;
    }
    // printf("malloc %p\n", read->rels);

    // Positions in the relationship array and read (0-indexed).
    size_t rels_pos = 0;
    size_t read_pos = 0;
    size_t next_rels_pos;
    size_t next_read_pos;
    // 5' and 3' ends of the read, in the read coordinates (1-indexed).
    size_t read_end5 = 1;
    size_t read_end3 = read->len;
    
    // Initialize the CIGAR operation.
    CigarOp cigar;
    init_cigar(&cigar, read->cigar);

    // Read the first operation from the pre-validated CIGAR string.
    if (get_next_cigar_op(&cigar))
        {assert(0);}
    
    // Read the entire CIGAR string one operation at a time.
    while (cigar.op != NULL)
    {
        assert(rels_pos <= read->num_rels);
        assert(read_pos <= read->len);
        // Decide what to do based on the current CIGAR operation.
        char op = *cigar.op;
        switch (op)
        {
            case CIG_ALIGN:
            case CIG_MATCH:
            case CIG_SUBST:
                // The read has matches or substitutions.
                next_rels_pos = rels_pos + cigar.len;
                assert(next_rels_pos <= read->num_rels);
                next_read_pos = read_pos + cigar.len;
                assert(next_read_pos <= read->len);
                while (rels_pos < next_rels_pos)
                {
                    read->rels[rels_pos] = encode_relate(rels_seq[rels_pos],
                                                         read->seq[read_pos],
                                                         read->quals[read_pos],
                                                         min_qual);
                    rels_pos++;
                    read_pos++;
                }
                break;
            
            case CIG_DELET:
                // Bases in the reference are deleted from the read.
                next_rels_pos = rels_pos + cigar.len;
                assert(next_rels_pos <= read->num_rels);
                // A deletion cannot occur first.
                if (rels_pos == 0)
                {
                    PyErr_SetString(RelateError,
                                    "A deletion was the first relationship");
                    return -1;
                }
                // If rels_pos > 0, then the only way for read_pos == 0
                // would be for the previous operation to have also been
                // a deletion, so an error about identical consecutive
                // CIGAR operations would have been raised already.
                assert(read_pos > 0);
                // A deletion cannot occur last.
                if (next_rels_pos == read->num_rels)
                {
                    PyErr_SetString(RelateError,
                                    "A deletion was the last relationship");
                    return -1;
                }
                // If next_rels_pos < read->num_rels, then the only way
                // for read_pos == read->len would be for all future
                // operations also to be deletions, so an error about
                // identical consecutive CIGAR operations would have
                // been raised already.
                assert(read_pos < read->len);
                // Record the deletion(s).
                while (rels_pos < next_rels_pos)
                {
                    // Mark the deletion in the array of relationships.
                    read->rels[rels_pos] = DELET;
                    // Add a deletion to the record of indels.
                    if (add_indel(&(read->pods), DELETION, rels_pos, read_pos))
                        {return -1;}
                    rels_pos++;
                }
                break;
            
            case CIG_INSRT:
                // Bases are inserted into the read.
                next_read_pos = read_pos + cigar.len;
                assert(next_read_pos <= read->len);
                if (rels_pos == 0)
                {
                    PyErr_SetString(RelateError,
                                    "An insertion was the first relationship");
                    return -1;
                }
                // If rels_pos > 0, then the only way for read_pos == 0
                // would be for the previous operation to have been a
                // deletion, so an error would have already been raised
                // by the deletion with read_pos == 0.
                assert(read_pos > 0);
                if (rels_pos == read->num_rels)
                {
                    PyErr_SetString(RelateError,
                                    "An insertion was the last relationship");
                    return -1;
                }
                // If rels_pos < read->num_rels, then the only way for
                // next_read_pos == read->len would be for all future
                // operations to be deletions, so an error would have
                // already been raised about insertions and deletions
                // being adjacent.
                assert(next_read_pos < read->len);
                // Record one insertion for each inserted base. Every
                // mutation needs a position in the relationship array.
                // But since an inserted base is (by definition) absent
                // from the reference, it does not correspond with one
                // position but lies between two positions. Either could
                // be used as the insertion's position; this algorithm
                // uses the 3' position. For example, if two bases are
                // inserted between reference positions 45 and 46, then
                // both inserted bases will be assigned to position 46.
                // Using the 3' position makes the math simpler than it
                // would be if using the 5' position because inserted
                // bases lie between the previous and subsequent CIGAR
                // operations, and ref_pos is the position immediately
                // 3' of the previous operation, so ref_pos is naturally
                // the position 3' of the insertion.
                while (read_pos < next_read_pos)
                {
                    // Add an insertion to the record of indels.
                    if (add_indel(&(read->pods), INSERTION, read_pos, rels_pos))
                        {return -1;}
                    read_pos++;
                }
                break;
            
            case CIG_SCLIP:
                // Bases were soft-clipped from the 5' or 3' end of the
                // read during alignment. Like insertions, they consume
                // the read but not the reference; however, as they are
                // not mutations, they do not need to be marked.
                next_read_pos = read_pos + cigar.len;
                assert(next_read_pos <= read->len);
                if (read_pos == 0)
                {
                    // This is the 5' soft clip.
                    // If read_pos == 0, then the only way that rels_pos
                    // could be > 0 is if there were deletions (and no
                    // other types of operations) before this one, but
                    // in that case the deletion would have raised an
                    // error already for being the first relationship.
                    assert(rels_pos == 0);
                    // If no previous consumed the reference or read
                    // (including soft clips), then this must be the
                    // only operation that can move read_end5.
                    assert(read_end5 == 1);
                    // Soft clip bases from the the 5' end of the read.
                    read_end5 += cigar.len;
                    assert(read_end5 <= read->len + 1);
                }
                else
                {
                    // This is the 3' soft clip.
                    // Check if any operations after this one consume
                    // the reference or read.
                    if (rels_pos < read->num_rels || next_read_pos < read->len)
                    {
                        PyErr_SetString(RelateError,
                                        "A soft clip occurred in the middle");
                        return -1;
                    }
                    // If no more operations consume the reference or
                    // read (including soft clips), then this must be
                    // the only operation that can move read_end3.
                    assert(read_end3 == read->len);
                    // Clip bases from the 3' end of the read.
                    assert(read_end3 == read_pos + cigar.len);
                    read_end3 = read_pos;  // equivalent to -= cigar.len
                }
                read_pos = next_read_pos;
                break;
            
            default:
                // The CIGAR string was already validated to contain no
                // other operations, so the default should never happen.
                assert(0);
        }

        // Read the next operation from the pre-validated CIGAR string.
        if (get_next_cigar_op(&cigar))
            {assert(0);}
        assert(cigar.op == NULL || *cigar.op != op);
    }

    // After reading all CIGAR operations, the positions in the read and
    // relationship array must have advanced to the last positions.
    assert(rels_pos == read->num_rels);
    assert(read_pos == read->len);

    // Add insertions to rels.
    const unsigned char ins_rel = get_ins_rel(insert3);
    for (size_t p = 0; p < read->pods.num_pods; p++)
    {
        IndelPod *pod = &(read->pods.pods[p]);
        if (pod->insert)
        {
            for (size_t i = 0; i < pod->num_indels; i++)
            {
                size_t ins_pos = get_lateral(pod->indels[i].lateral3, insert3);
                if (ins_pos < read->num_rels)
                {
                    if (read->rels[ins_pos] == MATCH)
                    {
                        // Insertions override matches.
                        read->rels[ins_pos] = ins_rel;
                    }
                    else
                    {
                        // Any other type of relationship can be marked
                        // on the same position as an insertion.
                        read->rels[ins_pos] |= ins_rel;
                    }
                }
            }
        }
    }
    
    if (ambindel)
    {
        // Find and mark ambiguous insertions and deletions.
        find_ambindels(read,
                       rels_seq,
                       read_end5,
                       read_end3,
                       min_qual,
                       insert3);
    }

    return 0;
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
                            const unsigned char *rels,
                            size_t read_pos,
                            size_t end5,
                            size_t end3)
{
    // Validate the arguments; positions are 1-indexed, and only end3
    // is allowed to be 0.
    assert(rels != NULL);
    assert(read_pos >= 1);
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
                              const SamRead *fwd_read,
                              const SamRead *rev_read)
{
    // Validate the arguments; positions are 1-indexed and the 5' ends
    // must be ≥ 1.
    assert(fwd_read->ref_end5 >= 1);
    assert(rev_read->ref_end5 >= 1);

    // Find the region where both reads 1 and 2 overlap.
    // All positions are 1-indexed.
    size_t both_end5, both_end3;

    // Find the region 5' of the overlap.
    if (fwd_read->ref_end5 > rev_read->ref_end5)
    {
        // The forward read begins after the reverse read.
        both_end5 = fwd_read->ref_end5;
        if (put_rels_in_dict(rels_dict,
                             rev_read->rels,
                             rev_read->pos,
                             rev_read->ref_end5,
                             min(rev_read->ref_end3, both_end5 - 1)))
            {return -1;}
    }
    else
    {
        // The forward read begins with or before the reverse read.
        both_end5 = rev_read->ref_end5;
        if (put_rels_in_dict(rels_dict,
                             fwd_read->rels,
                             fwd_read->pos,
                             fwd_read->ref_end5,
                             min(fwd_read->ref_end3, both_end5 - 1)))
            {return -1;}
    }

    // Find the region 3' of the overlap.
    if (fwd_read->ref_end3 > rev_read->ref_end3)
    {
        // The forward read ends after the reverse read.
        both_end3 = rev_read->ref_end3;
        if (put_rels_in_dict(rels_dict,
                             fwd_read->rels,
                             fwd_read->pos,
                             max(fwd_read->ref_end5, both_end3 + 1),
                             fwd_read->ref_end3))
            {return -1;}
    }
    else
    {
        // The forward read ends with or before the reverse read.
        both_end3 = fwd_read->ref_end3;
        if (put_rels_in_dict(rels_dict,
                             rev_read->rels,
                             rev_read->pos,
                             max(rev_read->ref_end5, both_end3 + 1),
                             rev_read->ref_end3))
            {return -1;}
    }

    // Fill relationships in the region of overlap.
    assert(both_end5 >= fwd_read->pos);
    assert(both_end5 >= rev_read->pos);
    unsigned char rel;
    for (size_t pos = both_end5; pos <= both_end3; pos++)
    {
        rel = (fwd_read->rels[pos - fwd_read->pos]
               &
               rev_read->rels[pos - rev_read->pos]);
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




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Python Interface                                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




/* Clean up by freeing all dynamically allocated memory for the reads
and decrementing the reference count for py_object (if given).
Do NOT pass SamReads that have not yet been initialized with init_read()
to this function, because it will call free_read() on them. */
static PyObject *cleanup(SamRead *read1,
                         SamRead *read2,
                         PyObject *py_object)
{
    // Free the dynamically allocated memory for the reads.
    free_read(read1);
    free_read(read2);

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
    Py_ssize_t line1_len;
    const char *line2;
    Py_ssize_t line2_len;
    const char *ref;
    const char *ref_seq;
    Py_ssize_t ref_len;
    unsigned long min_mapq;
    unsigned char min_qual;
    int insert3;
    int ambindel;
    int overhangs;
    unsigned long clip_end5;
    unsigned long clip_end3;
    if (!PyArg_ParseTuple(args,
                          "s#s#ss#kbpppkk",
                          &line1,
                          &line1_len,
                          &line2,
                          &line2_len,
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

    // The length arguments of lines 1/2 must equal the real lengths.
    assert(strlen(line1) == (size_t)line1_len);
    assert(strlen(line2) == (size_t)line2_len);

    // Declare two SamReads to hold the results and initialize them to
    // prevent any errors caused by uninitialized pointers.
    SamRead read1, read2;
    init_read(&read1);
    init_read(&read2);

    // Determine whether the read is paired-end and, if so, whether it
    // is properly paired.
    Py_ssize_t num_mates;
    int paired, proper;
    if (*line2)
    {
        // Paired-end
        assert(line2_len > 0);
        num_mates = 2;
        paired = 1;
        // The read is properly paired if lines 1 and 2 are different
        // because the way that SEISMIC-RNA indicates that a read is
        // paired-end but improperly paired (i.e. has one mate) is by
        // passing identical lines for reads 1 and 2.
        proper = (line1_len != line2_len) || (strcmp(line1, line2) != 0);
    }
    else
    {
        // Single-end
        assert(line2_len == 0);
        num_mates = 1;
        paired = 0;
        proper = 0;
    }
    
    // Initialize Python objects to return the results.

    PyObject *ends_rels_tuple = PyTuple_New(2);
    if (ends_rels_tuple == NULL)
        {return cleanup(&read1, &read2, ends_rels_tuple);}
    
    PyObject *ends_tuple = PyTuple_New(2);
    if (ends_tuple == NULL)
        {return cleanup(&read1, &read2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_rels_tuple, 0, ends_tuple);

    PyObject *end5s_list = PyList_New(num_mates);
    if (end5s_list == NULL)
        {return cleanup(&read1, &read2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_tuple, 0, end5s_list);

    PyObject *end3s_list = PyList_New(num_mates);
    if (end3s_list == NULL)
        {return cleanup(&read1, &read2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_tuple, 1, end3s_list);
    
    PyObject *rels_dict = PyDict_New();
    if (rels_dict == NULL)
        {return cleanup(&read1, &read2, ends_rels_tuple);}
    PyTuple_SET_ITEM(ends_rels_tuple, 1, rels_dict);

    // Calculate relationships for line 1.
    if (calc_rels_line(&read1,
                       line1,
                       (size_t)line1_len,
                       ref,
                       ref_seq,
                       (size_t)ref_len,
                       min_mapq,
                       paired,
                       proper,
                       min_qual,
                       insert3,
                       ambindel,
                       clip_end5,
                       clip_end3))
        {return cleanup(&read1, &read2, ends_rels_tuple);}

    if (paired)
    {
        // The read comprises forward (fwd) and reverse (rev) mates.
        SamRead *fwd_read, *rev_read;
        // Check if line 2 differs from line 1.
        if (proper)
        {
            // Calculate relationships for line 2.
            if (calc_rels_line(&read2,
                               line2,
                               (size_t)line2_len,
                               ref,
                               ref_seq,
                               (size_t)ref_len,
                               min_mapq,
                               paired,
                               proper,
                               min_qual,
                               insert3,
                               ambindel,
                               clip_end5,
                               clip_end3))
                {return cleanup(&read1, &read2, ends_rels_tuple);}

            // Check if reads 1 and 2 are paired properly.
            if (validate_pair(&read1, &read2))
                {return cleanup(&read1, &read2, ends_rels_tuple);}

            // Determine which read is forward and which is reverse.
            if (read2.reverse)
            {
                fwd_read = &read1;
                rev_read = &read2;
            }
            else
            {
                fwd_read = &read2;
                rev_read = &read1;
            }

            // Optionally remove overhangs.
            if (!overhangs)
            {
                // The 5' end of the reverse mate cannot extend past
                // the 5' end of the forward mate.
                if (rev_read->ref_end5 < fwd_read->ref_end5)
                {
                    rev_read->ref_end5 = fwd_read->ref_end5;
                }
                // The 3' end of the forward mate cannot extend past
                // the 3' end of the reverse mate.
                if (fwd_read->ref_end3 > rev_read->ref_end3)
                {
                    fwd_read->ref_end3 = rev_read->ref_end3;
                }
            }

            // Merge reads 1 and 2.
            if (put_2_rels_in_dict(rels_dict, fwd_read, rev_read))
                {return cleanup(&read1, &read2, ends_rels_tuple);}
        }
        else
        {
            // Lines 1 and 2 are identical: no need to process line 2.
            fwd_read = &read1;
            rev_read = &read1;
            if (put_rels_in_dict(rels_dict,
                                 read1.rels,
                                 read1.pos,
                                 read1.ref_end5,
                                 read1.ref_end3))
                {return cleanup(&read1, &read2, ends_rels_tuple);}
        }

        // Fill in the lists of 5' and 3' ends.
        if (put_end_in_list(end5s_list, 0, fwd_read->ref_end5))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
        if (put_end_in_list(end5s_list, 1, rev_read->ref_end5))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 0, fwd_read->ref_end3))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 1, rev_read->ref_end3))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
    }
    else
    {
        // Line 2 does not exist.
        if (put_rels_in_dict(rels_dict,
                             read1.rels,
                             read1.pos,
                             read1.ref_end5,
                             read1.ref_end3))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
        
        // Fill int the lists of 5' and 3' ends.
        if (put_end_in_list(end5s_list, 0, read1.ref_end5))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
        if (put_end_in_list(end3s_list, 0, read1.ref_end3))
            {return cleanup(&read1, &read2, ends_rels_tuple);}
    }

    // Free the memory of rels1 and rels2, but not of ends_rels_tuple
    // which will be returned.
    cleanup(&read1, &read2, NULL);

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
    NULL,          // module documentation
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


int main(void)
{
    // This module must be called via its Python API, not run directly.
    return 0;
}

/***********************************************************************
*                                                                      *
* © Copyright 2024, the Rouskin Lab.                                   *
*                                                                      *
* This file is part of SEISMIC-RNA.                                    *
*                                                                      *
* SEISMIC-RNA is free software; you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* SEISMIC-RNA is distributed in the hope that it will be useful, but   *
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- *
* ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     *
* Public License for more details.                                     *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  *
*                                                                      *
***********************************************************************/
