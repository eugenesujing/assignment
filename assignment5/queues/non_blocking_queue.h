#include "../common/allocator.h"

#define LFENCE asm volatile("lfence" : : : "memory")
#define SFENCE asm volatile("sfence" : : : "memory")



template <class T>
struct pointer_t{
  T* ptr;
  T* address(){
    T* temp = ptr;
    intptr_t p = ((intptr_t)temp <<16) >>16;
    return (T*)p;
  }
  uint16_t count(){
    T* temp = ptr;
    return (uintptr_t)temp>>48;
  }
};



template <class T>
struct Node
{
  T value;
  pointer_t<Node<T>> next;
};

template <class T>
class NonBlockingQueue
{
public:
  pointer_t<Node<T>> q_head;
  pointer_t<Node<T>> q_tail;
    CustomAllocator my_allocator_;
public:

    NonBlockingQueue() : my_allocator_()
    {
        std::cout << "Using NonBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next.ptr = NULL;
        q_head.ptr = q_tail.ptr = newNode;
    }

    Node<T>* compute_ptr(Node<T>* node, uint16_t count){
      Node<T>* temp = node;
      const uintptr_t mask = ~(0ULL);
      Node<T>* temp2 = (Node<T>*)(((uintptr_t)temp) | ((uint64_t)count <<48));
      return temp2;
    }

    void enqueue(T value)
    {
        // Use LFENCE and SFENCE as mentioned in pseudocode
        Node<T>* node = (Node<T>* )my_allocator_.newNode();
        node->value = value;
        node->next.ptr = NULL;
        SFENCE;
        pointer_t<Node<T>> tail;
        while(true) {
            tail = q_tail;
            LFENCE;
            pointer_t<Node<T>> next = tail.address()->next;
            LFENCE;
            std::cout<<"enque here 0"<<std::endl;
            if (tail.ptr == q_tail.ptr){
                if (next.address() == NULL) {
                //  std::cout<<"enque here 1"<<std::endl;
                  Node<T>* p = compute_ptr(node,next.count()+1);
                  pointer_t<Node<T>> to_be_swapped ={p};
                    if(CAS(&tail.address()->next, next, to_be_swapped)){
                      //std::cout<<"successfully enque "<<value<<std::endl;
                      break;
                    }

                }
                else{
                //  std::cout<<"enque here 2"<<std::endl;
                  Node<T>* p = compute_ptr(next.address(),tail.count()+1);
                  pointer_t<Node<T>> to_be_swapped ={p};
                    CAS(&q_tail, tail, to_be_swapped);	// ELABEL
                }
            }
        }
        SFENCE;
        Node<T>* p = compute_ptr(node,tail.count()+1);
        pointer_t<Node<T>> to_be_swapped ={p};
        CAS(&q_tail, tail, to_be_swapped);
    }

    bool dequeue(T *value)
    {
        // Use LFENCE and SFENCE as mentioned in pseudocode
        pointer_t<Node<T>> head;
        while(true){
              head = q_head;
             LFENCE;
             pointer_t<Node<T>> tail = q_tail;
             LFENCE;
             pointer_t<Node<T>> next = head.address()->next;
             LFENCE;
             if (head.ptr == q_head.ptr) {
               //std::cout<<"deque here 0"<<std::end;
                 if(head.address() == tail.address()) {
                   //std::cout<<"deque here 1"<<std::endl;
                     if(next.address() == NULL)
                             return false;
                    Node<T>* p = compute_ptr(next.address(),tail.count()+1);
                    pointer_t<Node<T>> to_be_swapped ={p};
                     CAS(&q_tail, tail, to_be_swapped);	//DLABEL
                 }
                 else {
                   //std::cout<<"deque here 1"<<std::endl;
                     *value = next.address()->value;
                     Node<T>* p = compute_ptr(next.address(),head.count()+1);
                     pointer_t<Node<T>> to_be_swapped ={p};
                     if(CAS(&q_head, head, to_be_swapped)){
                      // std::cout<<"successfully deque "<<*value<<std::endl;
                       break;
                     }

                 }
             }
         }
         my_allocator_.freeNode(head.address());
         return true;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }

};
