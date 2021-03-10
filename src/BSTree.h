/////////////////////////////////////////////////////////////////////
// Binary Search Tree.                                             //
// Implements the AVL tree algorithm.                              //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- History ---------------------------------------------------- //
// 04/02/2020 - Created - Brennan Young                            //
// 04/06/2020 - Modified - Brennan Young                           //
// - corrected error in copy.                                      //
// 04/10/2020 - Modified - Brennan Young                           //
// - implemented iterators.                                        //
// 06/02/2020 - Modified - Brennan Young                           //
// - defined const versions of iterator functions begin, last,     //
//   end, and at.                                                  //
// - ensured NULL nodes are kept as NULL nodes when copying.       //
// - clear now correctly sets node pointers, including the root,   //
//   to NULL.                                                      //
// 06/04/2020 - Modified - Brennan Young                           //
// - iterator and const_iterator now have default constructors.    //
// 07/27/2020 - Modified - Brennan Young                           //
// - added N and size(), for tracking size of container.           //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STDLIB_BSTREE_20200402
#define YOUNG_STDLIB_BSTREE_20200402

#include <cstddef> // size_t

namespace bystd { // Brennan Young standard namespace

template <class T>
class BSTree {
private:
    // data node
    struct Node {
        int key;        // key to sort and find
        T data;         // stored element
        int n;          // height of tree at element
        Node * parent;  // back up the tree
        Node * left;    // less-than path
        Node * right;   // greater-than path
        
        Node(int k=0) : key(k), n(1), parent(NULL), left(NULL),
            right(NULL) {}
    }; // BSTree::Node
    
    size_t N;
    Node * root;
    Node * current;
    
public:
    // iterators
    class iterator {
    private:
        Node* ptr;
    public:
        iterator();
        iterator(Node*);
        iterator& operator++();
        iterator operator++(int);
        iterator& operator--();
        iterator operator--(int);
        T& operator*();
        T* operator->();
        int key() const;
        bool operator==(const iterator&);
        bool operator!=(const iterator&);
    }; //BSTree::iterator
    
    class const_iterator {
    private:
        Node* ptr;
    public:
        const_iterator();
        const_iterator(Node*);
        const_iterator& operator++();
        const_iterator operator++(int);
        const_iterator& operator--();
        const_iterator operator--(int);
        const T& operator*();
        const T* operator->();
        int key() const;
        bool operator==(const const_iterator&);
        bool operator!=(const const_iterator&);
    }; // BSTree::const_iterator
    
private:
    // helpers
    Node* copy(Node*, const Node*);
    Node* clear(Node*);
    int height(const Node*) const;
    int balance(const Node*) const;
    Node* leftRotate(Node*);
    Node* rightRotate(Node*);
    Node* insert(Node*, int, const T&);
    Node* remove(Node*, int);
    void printKeys(const Node*) const;
    
public:
    // constructors / destructor
    BSTree();
    BSTree(const BSTree<T>&);
    ~BSTree();
    
    // operators
    BSTree<T> & operator=(const BSTree<T>&);
    
    // getters / setters
    size_t size() const;
    int height() const;
    iterator find(int);
    const_iterator find(int) const;
    int getKey(const iterator&) const;
    void insert(int, const T&);
    void remove(int);
    void clear();
    void printKeys() const;
    
    // iteration
    iterator begin();
    iterator last();
    iterator end();
    iterator at(int);
    const_iterator begin() const;
    const_iterator last() const;
    const_iterator end() const;
    const_iterator at(int) const;
};

// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor.
template<class T>
BSTree<T>::BSTree () : N(0), root(NULL), current(NULL) {}

// Copy constructor.
template<class T>
BSTree<T>::BSTree ( const BSTree<T>& bst ) :
    N(bst.N), root(NULL), current(NULL)
{
    root = copy(root, bst.root);
}

// Destructor.
template<class T>
BSTree<T>::~BSTree () { clear(root); }

// HELPERS //////////////////////////////////////////////////////////

// Copy all elements of subtree begining with node.
template<class T>
typename BSTree<T>::Node* BSTree<T>::copy ( BSTree<T>::Node* node,
    const BSTree<T>::Node* other )
{
    if ( other == NULL ) return NULL;
    
    // this node
    node = new Node (other->key);
    node->n = other->n;
    node->data = other->data;
    
    // children
    node->left = copy(node->left, other->left);
    if ( node->left != NULL ) node->left->parent = node;
    node->right = copy(node->right, other->right);
    if ( node->right != NULL ) node->right->parent = node;
    
    return node;
}

// Remove all elements of subtree beginning with node.
template<class T>
typename BSTree<T>::Node* BSTree<T>::clear ( BSTree<T>::Node* node )
{
    if ( node == NULL ) return NULL;
    node->left = clear(node->left);
    node->right = clear(node->right);
    delete node;
    node = NULL;
    --N;
    return NULL;
}

// Determine the height of the  subtree.
template<class T>
int BSTree<T>::height ( const BSTree<T>::Node* node ) const
{
    if ( node == NULL ) return 0;
    int L = node->left == NULL ? 0 : node->left->n;
    int R = node->right == NULL ? 0 : node->right->n;
    return 1 + (L > R ? L : R);
}

// Calculate the balance of the  subtree.
template<class T>
int BSTree<T>::balance ( const BSTree<T>::Node* node ) const
{
    if ( node == NULL ) return 0;
    int L = node->left == NULL ? 0 : node->left->n;
    int R = node->right == NULL ? 0 : node->right->n;
    return L - R;
}

// Left-rotation operation, to maintain tree balance.
template<class T>
typename BSTree<T>::Node* BSTree<T>::leftRotate (
    BSTree<T>::Node* node )
{
    Node* y = node->right;
    Node* z = y->left;
    
    y->left = node;
    node->right = z;
    
    if ( z != NULL ) z->parent = node;
    y->parent = node->parent;
    node->parent = y;
    
    node->n = height(node);
    y->n = height(y);
    
    return y; // new root
}

// Right-rotation operation, to maintain tree balance.
template<class T>
typename BSTree<T>::Node* BSTree<T>::rightRotate (
    BSTree<T>::Node* node )
{
    Node* y = node->left;
    Node* z = y->right;
    
    y->right = node;
    node->left = z;
    
    if ( z != NULL ) z->parent = node;
    y->parent = node->parent;
    node->parent = y;
    
    node->n = height(node);
    y->n = height(y);
    
    return y; // new root
}

// Add an element, if one with the same key does not already exist.
template<class T>
typename BSTree<T>::Node* BSTree<T>::insert (
    BSTree<T>::Node* node, int key, const T& x )
{
    // new node
    if ( node == NULL ) {
        node = new Node(key);
        ++N;
        node->data = x;
        return node;
    }
    
    // left branch
    if ( key < node->key ) {
        node->left = insert(node->left, key, x);
        node->left->parent = node;
    }
    
    // right branch
    else if ( node->key < key ) {
        node->right = insert(node->right, key, x);
        node->right->parent = node;
    }
    
    // no duplicates
    else return node;
    
    node->n = height(node);
    
    // re-balance
    int b = balance(node);
    
    if ( b > 1 && key < node->left->key ) // L-L
        return rightRotate(node);
    else if ( b < -1 && key > node->right->key ) // R-R
        return leftRotate(node);
    else if ( b > 1 && key > node->left->key ) { // L-R
        node->left = leftRotate(node->left);
        return rightRotate(node);
    }
    else if ( b < -1 && key < node->right->key ) { // R-L
        node->right = rightRotate(node->right);
        return leftRotate(node);
    }
    
    return node;
}

// Remove the element with key, if it exists.
template <class T>
typename BSTree<T>::Node* BSTree<T>::remove (
    BSTree<T>::Node* node, int key )
{
    if ( node == NULL ) return node;
    
    // left branch
    if ( key < node->key )
        node->left = remove(node->left, key);
    
    // right branch
    else if ( node->key < key )
        node->right = remove(node->right, key);
    
    // the node to remove
    else {
        //if ( current == node ) next();
        
        // no children
        if ( node->left == NULL && node->right == NULL ) {
            delete node;
            node = NULL;
            --N;
        }
        
        // one child (left)
        else if ( node->left != NULL && node->right == NULL ) {
            Node* t = node->left;
            node->left->parent = node->parent;
            delete node;
            node = t;
            --N;
        }
        
        // one child (right)
        else if ( node->left == NULL && node->right != NULL ) {
            Node* t = node->right;
            node->right->parent = node->parent;
            delete node;
            node = t;
            --N;
        }
        
        // two children
        else {
            // replace this node with smallest in right subtree
            Node* t = node->right;
            while ( t->left != NULL ) t = t->left;
            node->key = t->key;
            node->data = t->data;
            node->right = remove(node->right, t->key);
        }
    }
    
    if ( node == NULL ) return node;
    
    node->n = height(node);
    
    // re-balance
    int b = balance(node);
    
    if ( b > 1 && balance(node->left) >= 0 ) // L-L
        return rightRotate(node);
    else if ( b < -1 && balance(node->right) <= 0 ) // R-R
        return leftRotate(node);
    else if ( b > 1 && balance(node->left) < 0 ) { // L-R
        node->left = leftRotate(node->left);
        return rightRotate(node);
    }
    else if ( b < -1 && balance(node->right) > 0 ) { // R-L
        node->right = rightRotate(node->right);
        return leftRotate(node);
    }
    
    return node;
}

// Print the keys for each element in the tree (in order).
template<class T>
void BSTree<T>::printKeys ( const Node* node ) const
{
    if ( node == NULL ) return;
    printKeys(node->left);
    std::cout << node->key << " ";
    printKeys(node->right);
}

// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
template<class T>
BSTree<T> & BSTree<T>::operator= ( const BSTree<T> & bst )
{
    if ( &bst == this ) return *this; // copying self
    
    clear();
    N = bst.N;
    root = copy(root, bst.root);
    
    return *this;
}

// GETTERS / SETTERS ////////////////////////////////////////////////

// Return the number of items in the tree.
template<class T>
size_t BSTree<T>::size () const
{
    return N;
}

// Return the height of the tree.
template<class T>
int BSTree<T>::height () const
{
    return root == NULL ? 0 : root->n;
}

// Returns the item with key.
template<class T>
typename BSTree<T>::iterator BSTree<T>::find ( int key )
{
    Node* node = root;
    while ( node != NULL ) {
        if ( key < node->key ) node = node->left;
        else if ( node->key < key ) node = node->right;
        else return iterator(node);
    }
    return iterator(NULL);
}

template<class T>
typename BSTree<T>::const_iterator BSTree<T>::find ( int key ) const
{
    Node* node = root;
    while ( node != NULL ) {
        if ( key < node->key ) node = node->left;
        else if ( node->key < key ) node = node->right;
        else return const_iterator(node);
    }
    return const_iterator(NULL);
}

// Get the key of the element that the iterator is pointing at.
template<class T>
int BSTree<T>::getKey ( const BSTree<T>::iterator& it ) const
{
    return it.key();
}

// Add an element, if one with the same key does not already exist.
template<class T>
void BSTree<T>::insert ( int key, const T& x )
{
    root = insert(root, key, x);
}

// Remove an element, if one with the given key it exists.
template<class T>
void BSTree<T>::remove ( int key ) { root = remove(root, key); }

// Remove all elements.
template<class T>
void BSTree<T>::clear () {
    root = clear(root);
    current = NULL;
}

// Print the keys of all elements in the tree, in order.
template<class T>
void BSTree<T>::printKeys () const { printKeys(root); }

// ITERATION ////////////////////////////////////////////////////////

template <class T>
typename BSTree<T>::iterator BSTree<T>::begin ()
{
    Node* node = root;
    if ( node == NULL ) return iterator(node);
    while ( node->left != NULL ) node = node->left;
    return iterator(node);
}

template <class T>
typename BSTree<T>::iterator BSTree<T>::last ()
{
    Node* node = root;
    if ( node == NULL ) return iterator(node);
    while ( node->right != NULL ) node = node->right;
    return iterator(node);
}

template <class T>
typename BSTree<T>::iterator BSTree<T>::end ()
{
    return iterator(NULL);
}

template <class T>
typename BSTree<T>::iterator BSTree<T>::at ( int key )
{
    Node* node = root;
    while ( node != NULL ) {
        if ( key < node->key ) node = node->left;
        else if ( key > node->key ) node = node->right;
        else break;
    }
    return iterator(node);
}

template <class T>
typename BSTree<T>::const_iterator BSTree<T>::begin () const
{
    Node* node = root;
    if ( node == NULL ) return const_iterator(node);
    while ( node->left != NULL ) node = node->left;
    return const_iterator(node);
}

template <class T>
typename BSTree<T>::const_iterator BSTree<T>::last () const
{
    Node* node = root;
    if ( node == NULL ) return const_iterator(node);
    while ( node->right != NULL ) node = node->right;
    return const_iterator(node);
}

template <class T>
typename BSTree<T>::const_iterator BSTree<T>::end () const
{
    return const_iterator(NULL);
}

template <class T>
typename BSTree<T>::const_iterator BSTree<T>::at ( int key ) const
{
    Node* node = root;
    while ( node != NULL ) {
        if ( key < node->key ) node = node->left;
        else if ( key > node->key ) node = node->right;
        else break;
    }
    return const_iterator(node);
}


/////////////////////////////////////////////////////////////////////

// ITERATOR /////////////////////////////////////////////////////////

// Default constructor.
template<typename T>
BSTree<T>::iterator::iterator () : ptr(NULL) {}

// Constructor.
template<typename T>
BSTree<T>::iterator::iterator ( BSTree<T>::Node* it ) : ptr(it) {}

// Pre-increment.
template<typename T>
typename BSTree<T>::iterator& BSTree<T>::iterator::operator++ ()
{
    if ( ptr == NULL ) return *this;
    int key = ptr->key;
    
    // navigate to next root within which is the lowest-value node
    if ( ptr->right != NULL ) ptr = ptr->right;
    else while ( ptr != NULL && ptr->key <= key )
        ptr = ptr->parent;
    if ( ptr == NULL ) return *this;
    
    // get minimum-value node greater than key from root
    while ( ptr->left != NULL && ptr->left->key > key )
        ptr = ptr->left;
    
    return *this;
}

// Post-increment.
template<typename T>
typename BSTree<T>::iterator BSTree<T>::iterator::operator++ (
    int dummy )
{
    BSTree<T>::iterator it = *this;
    ++(*this);
    return it;
}

// Pre-decrement
template<typename T>
typename BSTree<T>::iterator& BSTree<T>::iterator::operator-- ()
{
    if ( ptr == NULL ) return *this;
    int key = ptr->key;
    
    // navigate to next root within which is the highest-value node
    if ( ptr->left != NULL ) ptr = ptr->left;
    else while ( ptr != NULL && ptr->key >= key )
        ptr = ptr->parent;
    if ( ptr == NULL ) return *this;
    
    // get maximum value node less than key from root
    while ( ptr->right != NULL && ptr->right->key < key )
        ptr = ptr->right;
    
    return *this;
}

// Post-decrement.
template<typename T>
typename BSTree<T>::iterator BSTree<T>::iterator::operator-- (
    int dummy )
{
    BSTree<T>::iterator it = *this;
    --(*this);
    return it;
}

// Data access.
template<typename T>
T& BSTree<T>::iterator::operator* ()
{
    return ptr->data;
}

template<typename T>
int BSTree<T>::iterator::key () const
{
    return ptr->key;
}

template<typename T>
T* BSTree<T>::iterator::operator-> ()
{
    return &ptr->data;
}

// Equivalence test.
template<typename T>
bool BSTree<T>::iterator::operator== ( 
    const BSTree<T>::iterator& it )
{
    return ptr == it.ptr;
}

template<typename T>
bool BSTree<T>::iterator::operator!= (
    const BSTree<T>::iterator& it )
{
    return ptr != it.ptr;
}


/////////////////////////////////////////////////////////////////////

// CONST ITERATOR ///////////////////////////////////////////////////

// Default constructor.
template<typename T>
BSTree<T>::const_iterator::const_iterator () : ptr(NULL) {}

// Constructor.
template<typename T>
BSTree<T>::const_iterator::const_iterator ( BSTree<T>::Node* it ) :
    ptr(it) {}

// Pre-increment.
template<typename T>
typename BSTree<T>::const_iterator&
    BSTree<T>::const_iterator::operator++ ()
{
    if ( ptr == NULL ) return *this;
    int key = ptr->key;
    
    // navigate to next root within which is the lowest-value node
    if ( ptr->right != NULL ) ptr = ptr->right;
    else while ( ptr != NULL && ptr->key <= key )
        ptr = ptr->parent;
    if ( ptr == NULL ) return *this;
    
    // get minimum-value node greater than key from root
    while ( ptr->left != NULL && ptr->left->key > key )
        ptr = ptr->left;
    
    return *this;
}

// Post-increment.
template<typename T>
typename BSTree<T>::const_iterator
    BSTree<T>::const_iterator::operator++ ( int dummy )
{
    BSTree<T>::const_iterator it = *this;
    ++(*this);
    return it;
}

// Pre-decrement
template<typename T>
typename BSTree<T>::const_iterator&
    BSTree<T>::const_iterator::operator-- ()
{
    if ( ptr == NULL ) return *this;
    int key = ptr->key;
    
    // navigate to next root within which is the highest-value node
    if ( ptr->left != NULL ) ptr = ptr->left;
    else while ( ptr != NULL && ptr->key >= key )
        ptr = ptr->parent;
    if ( ptr == NULL ) return *this;
    
    // get maximum value node less than key from root
    while ( ptr->right != NULL && ptr->right->key < key )
        ptr = ptr->right;
    
    return *this;
}

// Post-decrement.
template<typename T>
typename BSTree<T>::const_iterator
    BSTree<T>::const_iterator::operator-- ( int dummy )
{
    BSTree<T>::const_iterator it = *this;
    --(*this);
    return it;
}

// Data access.
template<typename T>
const T& BSTree<T>::const_iterator::operator* ()
{
    return ptr->data;
}

template<typename T>
int BSTree<T>::const_iterator::key () const
{
    return ptr->key;
}

template<typename T>
const T* BSTree<T>::const_iterator::operator-> ()
{
    return &ptr->data;
}

// Equivalence test.
template<typename T>
bool BSTree<T>::const_iterator::operator== ( 
    const BSTree<T>::const_iterator& it )
{
    return ptr == it.ptr;
}

template<typename T>
bool BSTree<T>::const_iterator::operator!= (
    const BSTree<T>::const_iterator& it )
{
    return ptr != it.ptr;
}

//

} // bystd

#endif // YOUNG_STDLIB_BSTREE_20200402