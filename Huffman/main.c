//
//  main.c
//  Huffman
//
//  Created by Locke Phan on 4/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#include <CoreFoundation/CoreFoundation.h>
#include <dispatch/dispatch.h>

#define DEBUG_ENABLED 0

void * allocate ( size_t size ) {
    void *space = malloc(size);
    if ( space == NULL ) {
        perror("Can't malloc");
        exit(EXIT_FAILURE);
    }
    return space;
}

void release ( void * data ) {
    free(data);
}

static char *code_word[256];
static size_t code_word_length[256];

struct Heap_Node {
    size_t              count;  
    int                 character_number;
    struct Heap_Node    *left;
    struct Heap_Node    *right;
};
typedef struct Heap_Node Heap_Node;

struct Stack_Node {
    char                value;
    struct Stack_Node   *previous;
    struct Stack_Node   *next;
};
typedef struct Stack_Node Stack_Node;

struct custom_stack_t {
    struct Stack_Node   *head;
    struct Stack_Node   *current;
    size_t              count;
};


/*************************************************************
 Push 'c' onto the stack
 
 Parameters:
 stack: pointer to the stack
 c: character to push on the stack
 *************************************************************/
void stack_push ( struct custom_stack_t *stack , char c ) {
    Stack_Node *newNode = (Stack_Node *) allocate(sizeof(Stack_Node));
    
    newNode->previous   = NULL;
    newNode->next       = NULL;
    newNode->value      = c;
    
    if (stack->head == NULL) {
        stack->head     = newNode;
        stack->current  = newNode;
    } else {
        stack->current->next    = newNode;
        newNode->previous       = stack->current;
        stack->current          = newNode;
    }
    stack->count++;
}

/*************************************************************
 Pop a node off the stack
 
 Parameters:
 stack: pointer to the stack
 *************************************************************/
void stack_pop ( struct custom_stack_t *stack ) {
    if (stack->current == NULL) return;    
    
    Stack_Node *popped_Node = stack->current;
    
    if (stack->current == stack->head) {
        stack->current  = NULL;
        stack->head     = NULL;
    } else {
        stack->current          = stack->current->previous;
        stack->current->next    = NULL;
    }
    release(popped_Node);
    popped_Node = NULL; // just to be safe
    stack->count--;
}

/************************************************************
 Returns a snapshot of the stack from top to bottom
 
 Parameters:
 stack: pointer to the stack
*************************************************************/
char * stack_print ( const struct custom_stack_t *stack ) {
    char *string = (char *) allocate(stack->count + 1);
    
    Stack_Node *p = stack->head;
    for(int i = 0 ; p != NULL ; p = p->next, i++)
        string[i] = p->value;
    
    string[stack->count] = '\0';
    return string;
}

/*************************************************************
 Pre-order traverse the Huffman Code tree to collect the code
 words
 
 Parameters:
 stack: pointer to a stack used to keep track of our code word
 root: pointer to the root node of a Huffman Code tree
**************************************************************/
void walk ( struct custom_stack_t *stack, const Heap_Node *root ) {
    if ( root == NULL ) return;
    if ( root->character_number >= 0 ) {
        char *myString = stack_print(stack);
        code_word[root->character_number] = myString;
        size_t length = strlen(myString);
        assert(1 <= length && length <= 8 && "Must be between 1 and 8 bits");
        code_word_length[root->character_number] = length;
    } else {
        stack_push(stack, '0');
        walk(stack, root->left);
        stack_pop(stack);
        
        stack_push(stack, '1');
        walk(stack, root->right);
        stack_pop(stack);
    }
}

/*************************************************************
 Compare predicate used by CoreFoundation Binary Heap
**************************************************************/
CFComparisonResult compare ( const void *ptr1, const void *ptr2, void *info ) {
    Heap_Node *node1 = (Heap_Node *) ptr1;
    Heap_Node *node2 = (Heap_Node *) ptr2;
    
    if ( node1->count > node2->count )
        return kCFCompareGreaterThan;
    else if ( node1->count < node2->count )
        return kCFCompareLessThan;
    else 
        return kCFCompareEqualTo;
}

/*************************************************************
 Creates a Huffman Code tree based on 'counts'
 
 Parameters:
 counts: pointer to an array of ints that describes the
 distribution of 256 bytes on a particular data set
**************************************************************/
Heap_Node * create_huffman_code_tree ( const int *counts ) {
    // Initializes binary heap to be used as a priority queue
    CFBinaryHeapCallBacks myCallBacks;
    myCallBacks.compare = compare;
    myCallBacks.retain = NULL;
    myCallBacks.release = NULL;
    CFBinaryHeapRef heap = CFBinaryHeapCreate(kCFAllocatorDefault, 0, &myCallBacks, NULL);
    
    // Loop over the array, create nodes, and
    // add them to the priority queue
    size_t total = 0;
    for (int i = 0 ; i < 256 ; i++) {
        if (counts[i] > 0) {
            total = total + counts[i];
            Heap_Node *heap_node        = (Heap_Node *) allocate(sizeof(Heap_Node));
            heap_node->count            = counts[i];
            heap_node->character_number = i;
            heap_node->left             = NULL;
            heap_node->right            = NULL;
            
            CFBinaryHeapAddValue(heap, heap_node);
        }
    }
    
    // Run Huffman Code algorithm to construct a prefix-free code tree
    CFIndex queue_size = CFBinaryHeapGetCount(heap) - 1;
    for ( int i = 0 ; i < queue_size; i++ ) {
        Heap_Node *node = (Heap_Node *) allocate(sizeof(Heap_Node));
        
        Heap_Node *x = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
        CFBinaryHeapRemoveMinimumValue(heap);
        Heap_Node *y = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
        CFBinaryHeapRemoveMinimumValue(heap);
        
        node->left = x;
        node->right = y;
        node->character_number = -1;
        node->count = x->count + y->count;
        
        CFBinaryHeapAddValue(heap, node);
    }
    
    // Extract the huffman code tree from the heap
    Heap_Node *node = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
    CFRelease(heap);
    
    return node;
}

/*************************************************************
 Recursively free a Huffman tree
 
 Parameters:
 node: pointer to the root node of the tree
*************************************************************/
void free_huffman_code_tree ( Heap_Node *node ) {
    Heap_Node *left = node->left;
    Heap_Node *right = node->right;
    release(node);
    
    if ( left != NULL )
        free_huffman_code_tree(left);
    if ( right != NULL )
        free_huffman_code_tree(right);
}


/*************************************************************
 Encodes 'data' and put the bytes into 'encodedData'
 
 Parameters:
 data: original uncompressed data
 encodedData: compressed data
 org_data_size: size of the original data (in bytes)
**************************************************************/
void encode (const unsigned char *data, unsigned char *encodedData, size_t org_data_size) {
    unsigned char byte = 0;
    size_t result_position = 0;
    size_t bit_position = 7;
    
    for ( size_t i = 0 ; i < org_data_size ; i++ ) {
        unsigned char c = data[i];
        
        char *code = code_word[c];
        
        size_t code_length = code_word_length[c];
        for ( int j = 0 ; j < code_length ; j++ ) {
            if ( code[j] == '1' ) 
                byte |= ( 1 << bit_position );
            if ( bit_position == 0 ) {
                encodedData[result_position] = byte;
                result_position++;
                bit_position = 7;
                byte = 0;
                continue;
            }
            bit_position--;
        }
    }
    if ( byte != 0 )
        encodedData[result_position] = byte;
    
}

/**********************************************************
 Decodes data using the tree that 'root' points to.
 Decoded data will be written into array 'decoded_data'
 
 Parameters:
 root: pointer to a Huffman Code tree
 decoded_data: pointer to an array of bytes
 encoded_data: pointer to an array of encoded bytes
 maxbits: number of bits encoded
**********************************************************/
void decode (Heap_Node *root, unsigned char *decoded_data, const unsigned char *encoded_data, size_t maxbits ) {    
    unsigned char byte;
    Heap_Node *current = root;
    size_t totalbits = 0;
    size_t decode_data_index = 0;
    size_t num_encoded_bytes = maxbits / 8;
    
    for ( size_t i = 0 ; i <= num_encoded_bytes ; i++ ) {
        byte = encoded_data[i];
        
        for ( int j = 7 ; j >= 0 && totalbits < maxbits ; j--, totalbits++ ) {
            if ( byte & ( 1 << j ))
                current = current->right;
            else
                current = current->left;
            if ( current->character_number >= 0 ) {
                decoded_data[decode_data_index] = current->character_number;
                current = root;
                decode_data_index++;
            }
        }
    }
}

/*****************************************************
 Returns the file size that 'path' points to in bytes
 
 Parameters:
 path: path to the file
******************************************************/
ssize_t get_file_size ( const char *path ) {
    struct stat buffer;
    int error;
    if ( (error = stat(path, &buffer)) < 0) {
        perror("Function: get_file_size.");
        exit(EXIT_FAILURE);
    }
    
    return buffer.st_size;
}

/*******************************************************
 Reads 'filesize' number of bytes using 'path' to 'buffer'
 Returns number of bytes read
 
 Parameters:
 path: path to the file
 filesize: size of the file
 buffer: a pointer to the buffer (expect to have 'filesize'
 number of bytes)
********************************************************/
ssize_t read_file_to_buffer ( const char *path, void *buffer, size_t filesize ) {
    int fd;
    if ((fd = open(path, O_RDONLY)) < 0) {
        perror("Can't open file");
        exit(EXIT_FAILURE);
    }
    
    ssize_t bytes_read;
    if ((bytes_read = read(fd, buffer, filesize)) < 0) {
        perror("Can't read read");
        exit(EXIT_FAILURE);
    }
    close(fd);
    return bytes_read;
}

/*******************************************************************
  Counts all the occurences of each char (or byte) using GCD.
 
 Parameters: 
    array : pointer to the array that will be updated with counts
            array is expected to have 1024 bytes in size
    data  : pointer to the data that will be counted
    data_size : size of the data in bytes
********************************************************************/
void count_dispatch ( int* array , const unsigned char *data , size_t data_size ) {
    // Zero out every single byte
    memset(array, 0, sizeof(int)*256); 

    // Determine the number of blocks using the size of data
    size_t NUM_BLOCKS = 1;
    for ( size_t i = 2 ; i <= 10 ; i++ ) { 
        if ( data_size % i == 0 ) {
            NUM_BLOCKS = i;
            break;
        }
    }

    size_t NUM_BYTE_PER_BLOCK = data_size / NUM_BLOCKS;
    
    // Create serial queue to access array
    dispatch_queue_t counts_Squeue = dispatch_queue_create("array", NULL);
    
    // Get global queue
    dispatch_queue_t global = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_group_t group = dispatch_group_create();
    
    clock_t time1 = clock();
    
    size_t upper_bound, lower_bound;
    lower_bound = 0;
    printf("Number of bytes per block: %ld. Block num: %ld.\n", NUM_BYTE_PER_BLOCK, NUM_BLOCKS);
    for ( size_t i = 1; i <= NUM_BLOCKS ; i++ ) {
        upper_bound = NUM_BYTE_PER_BLOCK * i;
        dispatch_group_async(group, global, ^{
            int *local_counts = (int *) allocate(sizeof(int)*256);
            memset(local_counts, 0, sizeof(int)*256);
            unsigned char c;
            for ( size_t j = lower_bound; j < upper_bound ; j++ ) {
                c = data[j];
                ++local_counts[c];
            }
            // Update local data to global count array
            dispatch_sync(counts_Squeue, ^{
                for ( int k = 0 ; k < 256 ; k++ )
                    array[k] += local_counts[k];
                release(local_counts);
            });
        });
        lower_bound = upper_bound;
    }
    
    dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
    
    clock_t time2 = clock();
    
    dispatch_release(group);
    dispatch_release(counts_Squeue);
    
    printf("Counting took %.4lf seconds using GCD.\n", (time2-time1)/(double)CLOCKS_PER_SEC);
}

/**************************************************************
 Counts all the occurences of each char (or byte) using a loop.
 
 Parameters: 
 array : pointer to the array that will be updated with counts
 array is expected to have 1024 bytes in size
 data  : pointer to the data that will be counted
 data_size : size of the data in bytes
**************************************************************/
void count_sequential ( int *array, const unsigned char *data, size_t data_size ) {
    memset(array, 0, sizeof(int)*256); // zero out every single byte
    
    // Counts the occurences of each byte using a simple loop
    clock_t time1 = clock();
    unsigned char c;
    for ( size_t i = 0; i < data_size ; i++ ) {
        c = data[i];
        ++array[c];
    }
    clock_t time2 = clock();
    
    printf("Counting took %.4lf seconds using traditional looping.\n", (time2-time1)/(double)CLOCKS_PER_SEC);
}

/**********************************************
 Uses counts to construct a Huffman Code tree.
 Collects and populate data into char *code_word[]
 And returns a pointer to the original code tree
 for decoding purposes.
 
 Parameters:
 counts: pointer to an array of counts
**********************************************/
Heap_Node* collect_code_words ( const int* counts ) {
    Heap_Node *node = create_huffman_code_tree(counts);
    
    // Initialize a stack to keep track of where we are in the tree
    struct custom_stack_t *stack = (struct custom_stack_t *) allocate(sizeof(struct custom_stack_t));
    stack->head = NULL;
    stack->current = NULL;
    stack->count = 0;
    
    // Traverse through the tree to collect our codewords
    walk(stack, node);
    
    return node;
}

int main (int argc, const char * argv[])
{
       char *filepath = "/Users/tienloc47/1meg.txt";
    //    char *filepath = "/Users/tienloc47/movie.avi";
    
    ssize_t filesize = get_file_size(filepath);
    printf("File size: %ld.\n", filesize);
    unsigned char *ORG_DATA = (unsigned char *) allocate(filesize);
    printf("Read %ld bytes into buffer.\n", read_file_to_buffer(filepath, ORG_DATA, filesize));
    
    // Initializes counts array to store the number of each character
    int *dispatch_counts = (int *) allocate(sizeof(int)*256);
    int *sequential_counts = (int *) allocate(sizeof(int)*256);
    
    // Counts the occurences of each byte using GCD
    count_dispatch(dispatch_counts, ORG_DATA, filesize);
    
    // Counts using traditional looping
    count_sequential(sequential_counts, ORG_DATA, filesize);
    
    // Compares gcd count with sequential count
    for ( int i = 0 ; i < 256 ; i++ ) 
        if ( sequential_counts[i] != dispatch_counts[i] )
            printf("%d: counts: %d, dispatchcounts: %d\n",i, sequential_counts[i], dispatch_counts[i]);

    // Collect code words
    Heap_Node *node = collect_code_words(dispatch_counts);
    
    // Loop over 256 bytes to count the length of each code word
    // and then sum them up to multiply with the count of each byte
    // to get the compressed size 
    size_t total_compressed_numbits = 0;
    for (int i = 0; i < 256 ; i++ ) {
        char *string = code_word[i];
        if (string == NULL) continue;
        size_t len = code_word_length[i];
        size_t size = len * sequential_counts[i];
        total_compressed_numbits += size;
    }
    
    printf("Size of file: %lu bytes.\n", node->count);
    printf("Size after compression: %lu bytes.\n", total_compressed_numbits/8);
    
    size_t compressed_byte_size = total_compressed_numbits / 8;
    
    unsigned char *compressed_data = (unsigned char *) allocate(compressed_byte_size);
    unsigned char *decompressed_data = (unsigned char *) allocate(filesize);
    memset(compressed_data, 0, compressed_byte_size);
    
#define SEQUENTIAL 1
#if SEQUENTIAL
    printf("\nStart sequential encode/decode...\n");
    clock_t time3 = clock();
    encode(ORG_DATA, compressed_data, filesize);
    decode(node, decompressed_data, compressed_data, total_compressed_numbits);
    clock_t time4 = clock();
    printf("    Sequential encode/decode took %.4lf seconds.\n", (time4-time3)/(double)CLOCKS_PER_SEC );
    
    printf("Comparing decoded data with original data...\n");
    int error = 0;
    for ( int i = 0; i < filesize ; i++ ) {
        if ( ORG_DATA[i] != decompressed_data[i] ) {
            error = 1;
            printf("Error! byte number: %d. ORG: %d , decoded: %d \n", i, ORG_DATA[i], decompressed_data[i]);
        }
    }
    if (!error)
        printf("Finished verifying.\n");
    else
        printf("Error occurred.\n");
#endif
    
#define PRODUCER_CONSUMER 1
#if PRODUCER_CONSUMER
    printf("\nStart producer-consumer encode/decode...\n");

    memset(compressed_data, 0, compressed_byte_size);
    memset(decompressed_data, 0, filesize);
    clock_t time1 = clock();

    dispatch_semaphore_t en_de_sem = dispatch_semaphore_create(0);
    dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
    dispatch_queue_t compressed_data_queue = dispatch_queue_create("compressed.data.access.queue", NULL);
    
    // Kick off the producer
    dispatch_async(global_queue, ^{
        __block unsigned char byte = 0;
        __block size_t result_position = 0;
        __block size_t bit_position = 7;
        
         void (^update_compress_data) (void) = ^{
             compressed_data[result_position] = byte; 
         };
        
        for ( size_t i = 0 ; i < filesize ; i++ ) {
            unsigned char c = ORG_DATA[i];  // Get a byte 
            char *code = code_word[c];      // Get the code word for that byte
            size_t code_length = code_word_length[c];
            for ( int j = 0 ; j < code_length ; j++ ) {
                if ( code[j] == '1' ) 
                    byte |= ( 1 << bit_position );
                if ( bit_position == 0 ) {
                    // Yo! I got 1 byte ready for ya
                    dispatch_sync(compressed_data_queue, update_compress_data);
                    // Increment the semaphore and signal
                    dispatch_semaphore_signal(en_de_sem);
                    result_position++;
                    bit_position = 7;
                    byte = 0;
                    continue;
                }
                bit_position--;
            }
        }
        // Write last byte to buffer
        if ( byte != 0 ) {
            dispatch_sync(compressed_data_queue, update_compress_data);
            dispatch_semaphore_signal(en_de_sem);
        }
    });
    
    // Wait for the consumer to finish
    dispatch_sync(global_queue, ^{
        __block unsigned char byte;
        Heap_Node *current = node;
        size_t totalbits = 0;
        size_t decode_data_index = 0;
        size_t num_encoded_bytes = total_compressed_numbits / 8;
        
        for ( size_t i = 0 ; i <= num_encoded_bytes && decode_data_index < filesize ; i++ ) {
            dispatch_semaphore_wait(en_de_sem, DISPATCH_TIME_FOREVER);
            dispatch_sync(compressed_data_queue, ^{ 
                byte = compressed_data[i]; 
            });
            
            for ( int j = 7 ; j >= 0 && totalbits < total_compressed_numbits ; j--, totalbits++ ) {
                if ( byte & ( 1 << j ))
                    current = current->right;
                else
                    current = current->left;
                if ( current->character_number >= 0 ) {
                    decompressed_data[decode_data_index] = current->character_number;
                    current = node;
                    decode_data_index++;
                }
            }
        }
    });
    
    clock_t time2 = clock();
    printf("    Producer-Consumer encode/decode took %.4lf seconds.\n", (time2-time1)/(double)CLOCKS_PER_SEC);
    
    printf("Comparing decoded data with original data...\n");
    error = 0;
    for ( size_t i = 0; i < filesize ; i++ ) {
        if ( ORG_DATA[i] != decompressed_data[i] ) {
            error = 1;
            printf("Error! byte number: %ld. ORG: %d , decoded: %d \n", i, ORG_DATA[i], decompressed_data[i]);
        }
    }
    if (!error)
        printf("Finished verifying.\n");
    else
        printf("Error!");
    
    dispatch_release(compressed_data_queue);
    dispatch_release(en_de_sem);
#endif    
    
    // Release encode array
    for ( int i = 0 ; i < 256 ; i++ )
        release(code_word[i]);
    release(sequential_counts);
    release(dispatch_counts);
    free_huffman_code_tree(node);
    
    printf("\nProgram finished executing.\n");
    return 0;
}

