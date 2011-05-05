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

void * allocate ( size_t size ) {
    void *space = malloc(size);
    if ( space == NULL ) {
        perror("Can't malloc");
        exit(EXIT_FAILURE);
    }
    return space;
}

char *code_word[256];

struct Heap_Node {
    size_t              count;  
    int                 character_number;
    struct Heap_Node    *left;
    struct Heap_Node    *right;
    struct Heap_Node    *parent;
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
    free(popped_Node);
    popped_Node = NULL; // just to be safe
    stack->count--;
}

char * stack_print ( struct custom_stack_t *stack ) {
    char *string = (char *) allocate(stack->count + 1);
    
    Stack_Node *p = stack->head;
    for(int i = 0 ; p != NULL ; p = p->next, i++)
        string[i] = p->value;
    
    string[stack->count] = '\0';
    return string;
}

void walk ( struct custom_stack_t *stack, Heap_Node *root ) {
    if ( root == NULL ) return;
    if ( root->character_number >= 0 ) {
        char *myString = stack_print(stack);
        code_word[root->character_number] = myString;
    } else {
        stack_push(stack, '0');
        walk(stack, root->left);
        stack_pop(stack);
        
        stack_push(stack, '1');
        walk(stack, root->right);
        stack_pop(stack);
    }
}


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


void free_huffman_code_tree ( Heap_Node *node ) {
    Heap_Node *left = node->left;
    Heap_Node *right = node->right;
    free(node);
    
    if ( left != NULL )
        free_huffman_code_tree(left);
    if ( right != NULL )
        free_huffman_code_tree(right);
}

Heap_Node * create_huffman_code_tree ( int *counts ) {
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
            heap_node->parent           = NULL;
            
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
        x->parent = node;
        node->right = y;
        y->parent = node;
        node->character_number = -1;
        node->count = x->count + y->count;
        
        CFBinaryHeapAddValue(heap, node);
    }
    
    // Extract the huffman code tree from the heap
    Heap_Node *node = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
    CFRelease(heap);
    
    return node;
}

void encode (unsigned char *data, unsigned char *encodedData, size_t ord_data_size) {
    unsigned char byte = 0;
    size_t result_position = 0;
    size_t bit_position = 7;
    
    for ( size_t i = 0 ; i < ord_data_size ; i++ ) {
        unsigned char c = data[i];
        
        char *code = code_word[c];
        
        size_t code_length = strlen(code);
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

void decode (Heap_Node *root, unsigned char *decoded_data, unsigned char *encoded_data, size_t maxbits ) {    
    unsigned char byte;
    Heap_Node *current = root;
    size_t totalbits = 0;
    size_t decode_data_index = 0;
    size_t num_encoded_bytes = maxbits / 8;
    
    for ( int i = 0 ; i <= num_encoded_bytes ; i++ ) {
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

ssize_t get_file_size ( char *path ) {
    struct stat buffer;
    int error;
    if ( (error = stat(path, &buffer)) < 0) {
        perror("Function: get_file_size.");
        exit(EXIT_FAILURE);
    }
    
    return buffer.st_size;
}

ssize_t read_file_to_buffer ( char *path, void *buffer, size_t filesize ) {
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
void count_dispatch ( int* array , unsigned char *data , size_t data_size ) {
    // Zero out every single byte
    memset(array, 0, sizeof(int)*256); 

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
        //printf("lower: %ld, upper: %ld.\n", lower_bound, upper_bound);
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
                free(local_counts);
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

int main (int argc, const char * argv[])
{
    //   char *filepath = "/Users/tienloc47/test.txt";
       char *filepath = "/Users/tienloc47/7meg.txt";
    //    char *filepath = "/Users/tienloc47/sample_large_text_file.txt";
    //    char *filepath = "/Users/tienloc47/movie.avi";
    //    char *filepath = "/Users/tienloc47/random_file";
#define DEBUGGING 0
    
    ssize_t filesize = get_file_size(filepath);
    printf("File size: %ld.\n", filesize);
    unsigned char *ORG_DATA = (unsigned char *) allocate(filesize);
    printf("Read %ld bytes into buffer.\n", read_file_to_buffer(filepath, ORG_DATA, filesize));
    
    // Initializes counts array to store the number of each character
    int *dispatch_counts = (int *) allocate(sizeof(int)*256);
    int *counts = (int *) allocate(sizeof(int)*256);
    memset(counts, 0, sizeof(int)*256); // zero out every single byte
    
    // Counts the occurences of each byte using GCD
    count_dispatch(dispatch_counts, ORG_DATA, filesize);
    
    // Counts the occurences of each byte using a simple loop
    clock_t time1 = clock();
    unsigned char c;
    for ( int i = 0; i < filesize ; i++ ) {
        c = ORG_DATA[i];
        ++counts[c];
    }
    clock_t time2 = clock();
    
    printf("Counting took %.4lf seconds using traditional looping.\n", (time2-time1)/(double)CLOCKS_PER_SEC);

    
    for ( int i = 0 ; i < 256 ; i++ ) {
        if ( counts[i] != dispatch_counts[i] )
            printf("%d: counts: %d, dispatchcounts: %d\n",i, counts[i], dispatch_counts[i]);
    }
    
    
    Heap_Node *node = create_huffman_code_tree(counts);
    
    // Initialize a stack to keep track of where we are in the tree
    struct custom_stack_t *stack = (struct custom_stack_t *) allocate(sizeof(struct custom_stack_t));
    stack->head = NULL;
    stack->current = NULL;
    stack->count = 0;
    
    // Traverse through the tree to collect our codewords
    walk(stack, node);
    
    // Loop over 256 bytes to count the length of each code word
    // and then sum them up to multiply with the count of each byte
    // to get the compressed size 
    size_t total_compressed_numbits = 0;
    for (int i = 0; i < 256 ; i++ ) {
        char *string = code_word[i];
        if (string == NULL) continue;
        size_t len = strlen(string);
        size_t size = len * counts[i];
        total_compressed_numbits += size;
    }
    
    printf("Size of file: %lu bytes.\n", node->count);
    printf("Size after compression: %lu bytes.\n", total_compressed_numbits/8);
    
    size_t compressed_byte_size = total_compressed_numbits / 8;
    
    unsigned char *compressed_data = (unsigned char *) allocate(compressed_byte_size);
    unsigned char *decoded_data = (unsigned char *) allocate(filesize);
    memset(compressed_data, 0, compressed_byte_size);
    
#define SEQUENTIAL 1
#if SEQUENTIAL
    clock_t time3 = clock();
    printf("Start encoding process...\n");
    time1 = clock();
    encode(ORG_DATA, compressed_data, filesize);
    time2 = clock();
    printf("Finished encoding. Took: %.4lf seconds.\n", (time2-time1)/(double)CLOCKS_PER_SEC);
    printf("Start decoding process...\n");
    time1 = clock();
    decode(node, decoded_data, compressed_data, total_compressed_numbits);
    time2 = clock();
    printf("Finished decoding. Took: %.4lf seconds.\n", (time2-time1)/(double)CLOCKS_PER_SEC);
    clock_t time4 = clock();
    printf("Sequential encode/decode took %.4lf seconds.\n", (time4-time3)/(double)CLOCKS_PER_SEC );
    
    printf("Comparing decoded data with original data...\n");
    for ( int i = 0; i < filesize ; i++ ) {
        if ( ORG_DATA[i] != decoded_data[i] )
            printf("Error! byte number: %d. ORG: %d , decoded: %d \n", i, ORG_DATA[i], decoded_data[i]);
    }
#endif
    
#define PRODUCER_CONSUMER 1
#if SEQUENTIAL && PRODUCER_CONSUMER
    memset(compressed_data, 0, compressed_byte_size);
    memset(decoded_data, 0, filesize);
#endif
#if PRODUCER_CONSUMER
    time1 = clock();
#if DEBUGGING
    __block size_t debug_sem = 0;
#endif
    dispatch_semaphore_t en_de_sem = dispatch_semaphore_create(0);
    dispatch_queue_t compressed_data_queue = dispatch_queue_create("compressed.data.access.queue", NULL);
    dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
        __block unsigned char byte = 0;
        __block size_t result_position = 0;
        __block size_t bit_position = 7;
        
         void (^update_compress_data) (void) = ^{
#if DEBUGGING
             printf("Writing 1 byte to buffer. Byte: %d.\n", byte);
#endif
             compressed_data[result_position] = byte; 
         };
        
        for ( size_t i = 0 ; i < filesize ; i++ ) {
            unsigned char c = ORG_DATA[i];
            char *code = code_word[c];
            size_t code_length = strlen(code);
            for ( int j = 0 ; j < code_length ; j++ ) {
                if ( code[j] == '1' ) 
                    byte |= ( 1 << bit_position );
                if ( bit_position == 0 ) {
                    // Yo! I got 1 byte ready for ya
                    dispatch_sync(compressed_data_queue, update_compress_data);
                    // Increment the semaphore and signal
                    dispatch_semaphore_signal(en_de_sem);
#if DEBUGGING
                    debug_sem++;
#endif
                    result_position++;
                    bit_position = 7;
                    byte = 0;
                    continue;
                }
                bit_position--;
            }
        }
        if ( byte != 0 ) {
            dispatch_sync(compressed_data_queue, update_compress_data);
            dispatch_semaphore_signal(en_de_sem);
#if DEBUGGING
            debug_sem++;
#endif
        }
#if DEBUGGING
        printf("I'm done encoding data! Sem: %ld.\n", debug_sem);
#endif
    });
    
    __block unsigned char byte;
    Heap_Node *current = node;
    size_t totalbits = 0;
    size_t decode_data_index = 0;
    size_t num_encoded_bytes = total_compressed_numbits / 8;
        
    for ( int i = 0 ; i <= num_encoded_bytes ; i++ ) {
        dispatch_semaphore_wait(en_de_sem, DISPATCH_TIME_FOREVER);
        dispatch_sync(compressed_data_queue, ^{ 
            byte = compressed_data[i]; 
#if DEBUGGING                
            printf("Read 1 byte from buffer. Read: %d.\n", byte);
#endif
        });
        
        for ( int j = 7 ; j >= 0 && totalbits < total_compressed_numbits ; j--, totalbits++ ) {
            if ( byte & ( 1 << j ))
                current = current->right;
            else
                current = current->left;
            if ( current->character_number >= 0 ) {
                decoded_data[decode_data_index] = current->character_number;
                current = node;
                decode_data_index++;
            }
        }
        if (decode_data_index == filesize)
            break;
    }
    dispatch_release(compressed_data_queue);
    dispatch_release(en_de_sem);
    time2 = clock();
    printf("Pro-Con encode/decode took %.4lf seconds.\n", (time2-time1)/(double)CLOCKS_PER_SEC);
    
    printf("Comparing decoded data with original data...\n");
    for ( int i = 0; i < filesize ; i++ ) {
        if ( ORG_DATA[i] != decoded_data[i] )
            printf("Error! byte number: %d. ORG: %d , decoded: %d \n", i, ORG_DATA[i], decoded_data[i]);
    }
#endif    
    
    // Release encode array
    for ( int i = 0 ; i < 256 ; i++ )
        free(code_word[i]);
    
    free(counts);    
    free(stack);
    free_huffman_code_tree(node);
    
    printf("Program finished executing.\n");
    return 0;
}

