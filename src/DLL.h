#ifndef YOUNG_STDLIB_LINKEDLIST_20180510
#define YOUNG_STDLIB_LINKEDLIST_20180510

/////////////////////////////////////////////////////////////////////
// 05/10/2018 - Created - Brennan Young                            //
// 05/24/2018 - Modified - Brennan Young                           //
// - removed get and set member functions--they are superfluous    //
//   with the [] operator.                                         //
// 08/02/2019 - Modified - Brennan Young                           //
// - operator [] const now returns const reference to item.        //
// - clear now starts at the root and moves through the list.      //
// 08/19/2019 - Modified - Brennan Young                           //
// - empty pointers are now NULL instead of 0.                     //
// 08/26/2019 - Modified - Brennan Young                           //
// - default constructor updated to use constructor : {} syntax.   //
// - minor format changes.                                         //
// 03/12/2020 - Modified - Brennan Young                           //
// - find, push_back, and insert now take a constant reference to  //
//   an object.                                                    //
// - pop renamed to pop_back, and pop_back doesn't return          //
//   anything.                                                     //
// - added push_front, pop_front.                                  //
// - pop_front and pop_back no longer return anything.             //
// - some minor improvements to style and efficiency.              //
// 03/31/2020 - Modified - Brennan Young                           //
// - added iteration capability: current, getCurrent, first, last, //
//   next, prev, and currentValid.                                 //
// 04/01/2020 - Modified - Brennan Young                           //
// - added insertBeforeCurrent and removeCurrent.                  //
// 04/07/2020 - Modified - Brennan Young                           //
// - modified iteration to return pointers.                        //
// - checked if current is NULL before attempting to move it       //
//   forward or backward.                                          //
// 04/09/2020 - Modified - Brennan Young                           //
// - corrected minor problems in push methods.                     //
// - initialized root and tail to NULL in copy constructor.        //
// 04/10/2020 - Modified - Brennan Young                           //
// - the current pointer has been adjusted to always point to the  //
///  next element if the current element is removed, but to        //
//   otherwise remain unchanged.                                   //
// - added private scan functions to standardize internal          //
//   iteration.                                                    //
// 05/13/2020 - Modified - Brennan Young                           //
// - replaced internal iteration with iterators.                   //
/////////////////////////////////////////////////////////////////////

#include <cstddef>

namespace bystd { // Brennan Young standard namespace

template <class T>
class DLL {
private:
    struct Elem {
        T value;
        Elem * next;
        Elem * back;
    };
    Elem * root;    // first element in list
    Elem * tail;    // last element in list
    size_t n; // number of elements
    
    // helpers
    Elem* scani(size_t) const;
    Elem* scan(const T&) const;
    Elem* scan(const T&, size_t*) const;
    
public:
    // iterators
    class iterator {
    private:
        Elem* ptr;
        friend class DLL<T>;
    public:
        iterator(Elem*);
        iterator& operator++();
        iterator operator++(int);
        iterator& operator--();
        iterator operator--(int);
        T& operator*();
        const T& operator*() const;
        T* operator->();
        const T* operator->() const;
        bool operator==(const iterator&) const;
        bool operator!=(const iterator&) const;
    }; //DLL::iterator
    
    class const_iterator {
    private:
        Elem* ptr;
        friend class DLL<T>;
    public:
        const_iterator(Elem*);
        const_iterator& operator++();
        const_iterator operator++(int);
        const_iterator& operator--();
        const_iterator operator--(int);
        const T& operator*() const;
        const T* operator->() const;
        bool operator==(const const_iterator&) const;
        bool operator!=(const const_iterator&) const;
    }; // DLL::const_iterator
    
public:
    // constructors, destructor
    DLL();
    DLL(const DLL<T>&);
    ~DLL();
    
    // operators
    DLL<T> & operator=(const DLL<T>&);
    T& operator[](size_t);
    const T& operator[](size_t) const;
    
    // getters
    size_t size() const;
    size_t find(const T&) const;
    size_t find(const iterator&) const;
    size_t find(const const_iterator&) const;
    
    // modifiers
    void push_front(const T&);
    void push_back(const T&);
    void insert(size_t, const T&);
    void insert(const iterator&, const T&);
    void insert(const const_iterator&, const T&);
    void pop_front();
    void pop_back();
    void remove(size_t);
    void remove(iterator&);
    void remove(const_iterator&);
    void clear();
    
    // iteration
    iterator begin();
    iterator last();
    iterator end();
    iterator at(size_t);
    const_iterator begin() const;
    const_iterator last() const;
    const_iterator end() const;
    const_iterator at(size_t) const;
};

// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor.
template<class T>
DLL<T>::DLL () : root(NULL), tail(NULL), n(0) {}

// Copy constructor.
template<class T>
DLL<T>::DLL ( const DLL<T> & L ) : root(NULL), tail(NULL)
{
    for ( Elem * Le = L.root; Le != NULL; Le = Le->next )
        this->push_back(Le->value);
}

// Destructor.
template<class T>
DLL<T>::~DLL () { clear(); }

// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
template<class T>
DLL<T> & DLL<T>::operator= ( const DLL<T> & L )
{
    if ( &L == this ) return *this;
    
    // clear the list
    this->clear();
    
    // add all the elements of L
    for ( Elem * Le = L.root; Le != NULL; Le = Le->next )
        this->push_back(Le->value);
    
    return *this;
}

// Element access.
template<class T>
T & DLL<T>::operator[] ( size_t i )
{
    Elem* e = scani(i);
    return e->value;
}

// Element access.
template<class T>
const T & DLL<T>::operator[] ( size_t i ) const
{
    Elem* e = scani(i);
    return e->value;
}

// HELPERS //////////////////////////////////////////////////////////

// Get the element at position i in the list.
template <class T>
typename DLL<T>::Elem* DLL<T>::scani ( size_t i ) const
{
    Elem* e = NULL;
    if ( i >= n ) return e;
    if ( i < 0.5 * n ) {
        e = root;
        for ( size_t j = 0; j < i; ++j ) e = e->next;
    }
    else {
        e = tail;
        for ( size_t j = n-1; j > i; --j ) e = e->back;
    }
    return e;
}

// Get the element with the given value.
template <class T>
typename DLL<T>::Elem* DLL<T>::scan ( const T& x ) const
{
    Elem* e = root;
    while ( e != NULL && !(e->value == x) ) e = e->next;
    return e;
}

// Get the element and index number of the element with the given
// value.
template <class T>
typename DLL<T>::Elem* DLL<T>::scan ( const T& x, size_t* i )
    const
{
    Elem* e = root;
    *i = 0;
    while ( e != NULL && !(e->value == x) ) {
        e = e->next;
        *i += 1;
    }
    return e;
}

// GETTERS / SETTERS ////////////////////////////////////////////////

template<class T>
size_t DLL<T>::size () const { return n; }

// Returns the index of the item, or the size of the list if the
// item is not found.
template<class T>
size_t DLL<T>::find ( const T & x ) const
{
    size_t i;
    Elem* e = scan(x, &i);
    return i;
}

// Returns the index corresponding to the iterator, or the size of
// the list if the iterator is not found.
template<class T>
size_t DLL<T>::find ( const DLL<T>::iterator& it ) const
{
    size_t i;
    Elem* e = root;
    while ( e != NULL && e != it.ptr ) e = e->next;
    return i;
}

// Returns the index corresponding to the iterator, or the size of
// the list if the iterator is not found.
template<class T>
size_t DLL<T>::find ( const DLL<T>::const_iterator& it ) const
{
    size_t i;
    Elem* e = root;
    while ( e != NULL && e != it.ptr ) e = e->next;
    return i;
}

// MODIFIERS ////////////////////////////////////////////////////////

// Add to front of list.
template<class T>
void DLL<T>::push_front ( const T & x )
{
    Elem * e = new Elem;
    e->value = x;
    e->next = root;
    e->back = NULL;
    if ( root != NULL ) root->back = e;
    root = e;
    if ( tail == NULL ) tail = e;
    ++n;
}

// Add to end of list.
template<class T>
void DLL<T>::push_back ( const T & x )
{
    Elem * e = new Elem;
    e->value = x;
    e->next = NULL;
    e->back = tail;
    if ( tail != NULL ) tail->next = e;
    tail = e;
    if ( root == NULL ) root = e;
    ++n;
}

// Remove from the front of the list.
template <class T>
void DLL<T>::pop_front ()
{
    if ( root == NULL ) return;
    Elem * e = root->next;
    delete root;
    root = e;
    if ( e == NULL ) tail = NULL;
    else e->back = NULL;
    --n;
}

// Remove from the end of the list.
template<class T>
void DLL<T>::pop_back ()
{
    if ( tail == NULL ) return;
    Elem * e = tail->back;
    delete tail;
    tail = e;
    if ( e == NULL ) root = NULL;
    else e->next = NULL;
    --n;
}

// Insert the item before the item at index i.
template<class T>
void DLL<T>::insert ( size_t i, const T& x ) {
    // insert at front
    if ( i == 0 ) {
        push_front(x);
        return;
    }
    
    // insert at back
    if ( i >= n ) {
        push_back(x);
        return;
    }
    
    // traverse to insert location
    Elem* t = scani(i);
    
    // insert new element
    Elem* e = new Elem;
    e->value = x;
    e->next = t;
    e->back = t->back;
    t->back->next = e;
    t->back = e;
    ++n;
}

// Insert the item at the iterator position.
template<typename T>
void DLL<T>::insert ( const DLL<T>::iterator& it, const T& x )
{
    if ( it == end() ) {
        push_back(x);
        return;
    }
    
    // insert new element
    Elem* e = new Elem;
    e->value = x;
    e->next = it.ptr;
    e->back = it.ptr->back;
    it.ptr->back->next = e;
    it.ptr->back = e;
    ++n;
}

// Insert the item at the iterator position.
template<typename T>
void DLL<T>::insert ( const DLL<T>::const_iterator& it, const T& x )
{
    if ( it == end() ) {
        push_back(x);
        return;
    }
    
    // insert new element
    Elem* e = new Elem;
    e->value = x;
    e->next = it.ptr;
    e->back = it.ptr->back;
    it.ptr->back->next = e;
    it.ptr->back = e;
    ++n;
}

// Remove the item at index i.
template<class T>
void DLL<T>::remove ( size_t i )
{
    if ( i >= n ) return;
    
    // traverse to remove location
    Elem* e = scani(i);
    
    // remove element
    if ( e->back != NULL ) e->back->next = e->next;
    else root = e->next;
    if ( e->next != NULL ) e->next->back = e->back;
    else tail = e->back;
    delete e;
    --n;
}

// Remove the item at the iterator.
template <class T>
void DLL<T>::remove ( DLL<T>::iterator& it )
{
    if ( it == end() ) return;
    
    // remove element
    if ( it.ptr->back != NULL ) it.ptr->back->next = it.ptr->next;
    else root = it.ptr->next;
    if ( it.ptr->next != NULL ) it.ptr->next->back = it.ptr->back;
    else tail = it.ptr->back;
    it.ptr = it.ptr->next;
    --n;
}

// Remove the item at the iterator.
template <class T>
void DLL<T>::remove ( DLL<T>::const_iterator& it )
{
    if ( it == end() ) return;
    
    // remove element
    if ( it.ptr->back != NULL ) it.ptr->back->next = it.ptr->next;
    else root = it.ptr->next;
    if ( it.ptr->next != NULL ) it.ptr->next->back = it.ptr->back;
    else tail = it.ptr->back;
    it.ptr = it.ptr->next;
    --n;
}

// Remove all items from the list.
template<class T>
void DLL<T>::clear ()
{
    Elem * e = root;
    Elem * ne;
    while ( e ) {
        ne = e->next;
        delete e;
        e = ne;
    }
    root = tail = NULL;
    n = 0;
}

// ITERATION ////////////////////////////////////////////////////////

template <class T>
typename DLL<T>::iterator DLL<T>::begin () { return iterator(root); }

template <class T>
typename DLL<T>::iterator DLL<T>::last () { return iterator(tail); }

template <class T>
typename DLL<T>::iterator DLL<T>::end () { return iterator(NULL); }

template <class T>
typename DLL<T>::iterator DLL<T>::at ( size_t i )
{ return iterator(scani(i)); }

template <class T>
typename DLL<T>::const_iterator DLL<T>::begin () const
{ return const_iterator(root); }

template <class T>
typename DLL<T>::const_iterator DLL<T>::last () const
{ return const_iterator(tail); }

template <class T>
typename DLL<T>::const_iterator DLL<T>::end () const
{ return const_iterator(NULL); }

template <class T>
typename DLL<T>::const_iterator DLL<T>::at ( size_t i ) const
{ return const_iterator(scani(i)); }


/////////////////////////////////////////////////////////////////////

// ITERATOR /////////////////////////////////////////////////////////

// Constructor.
template<typename T>
DLL<T>::iterator::iterator ( DLL<T>::Elem* it ) : ptr(it) {}

// Pre-increment.
template<typename T>
typename DLL<T>::iterator& DLL<T>::iterator::operator++ ()
{
    if ( ptr != NULL ) ptr = ptr->next;
    return *this;
}

// Post-increment.
template<typename T>
typename DLL<T>::iterator DLL<T>::iterator::operator++ (
    int dummy )
{
    DLL<T>::iterator it = *this;
    ++(*this);
    return it;
}

// Pre-decrement
template<typename T>
typename DLL<T>::iterator& DLL<T>::iterator::operator-- ()
{
    if ( ptr != NULL ) ptr = ptr->back;
    return *this;
}

// Post-decrement.
template<typename T>
typename DLL<T>::iterator DLL<T>::iterator::operator-- (
    int dummy )
{
    DLL<T>::iterator it = *this;
    --(*this);
    return it;
}

// Data access.
template<typename T>
T& DLL<T>::iterator::operator* () { return ptr->value; }

template<typename T>
const T& DLL<T>::iterator::operator* () const { return ptr->value; }

template<typename T>
T* DLL<T>::iterator::operator-> () { return &ptr->value; }

template<typename T>
const T* DLL<T>::iterator::operator-> () const
{ return &ptr->value; }

// Equivalence test.
template<typename T>
bool DLL<T>::iterator::operator== ( const DLL<T>::iterator& it )
    const
{ return ptr == it.ptr; }

template<typename T>
bool DLL<T>::iterator::operator!= ( const DLL<T>::iterator& it )
    const
{ return ptr != it.ptr; }


/////////////////////////////////////////////////////////////////////

// CONST ITERATOR ///////////////////////////////////////////////////

// Constructor.
template<typename T>
DLL<T>::const_iterator::const_iterator ( DLL<T>::Elem* it ) :
    ptr(it) {}

// Pre-increment.
template<typename T>
typename DLL<T>::const_iterator&
    DLL<T>::const_iterator::operator++ ()
{
    if ( ptr != NULL ) ptr = ptr->next;
    return *this;
}

// Post-increment.
template<typename T>
typename DLL<T>::const_iterator
    DLL<T>::const_iterator::operator++ ( int dummy )
{
    DLL<T>::const_iterator it = *this;
    ++(*this);
    return it;
}

// Pre-decrement
template<typename T>
typename DLL<T>::const_iterator&
    DLL<T>::const_iterator::operator-- ()
{
    if ( ptr != NULL ) ptr = ptr->back;
    return *this;
}

// Post-decrement.
template<typename T>
typename DLL<T>::const_iterator
    DLL<T>::const_iterator::operator-- ( int dummy )
{
    DLL<T>::const_iterator it = *this;
    --(*this);
    return it;
}

// Data access.
template<typename T>
const T& DLL<T>::const_iterator::operator* () const
{ return ptr->value; }

template<typename T>
const T* DLL<T>::const_iterator::operator-> () const
{ return &ptr->value; }

// Equivalence test.
template<typename T>
bool DLL<T>::const_iterator::operator== ( 
    const DLL<T>::const_iterator& it ) const
{ return ptr == it.ptr; }

template<typename T>
bool DLL<T>::const_iterator::operator!= (
    const DLL<T>::const_iterator& it ) const
{ return ptr != it.ptr; }

} // bystd

#endif // YOUNG_STDLIB_LINKEDLIST_20180510