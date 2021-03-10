/////////////////////////////////////////////////////////////////////
// ProgCounter                                                     //
// Object for tracking, managing, and reporting process progress.  //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 08/25/2017 - Brennan Young                                      //
// - created                                                       //
// modified up to 10/19/2017 - Brennan Young                       //
// 10/19/2017 - Brennan Young                                      //
// - added reset(ProgCounter).                                     //
// 11/15/2017 - Brennan Young                                      //
// - replaced the ProgCounter struct with the ProgCounter class.   //
// - overloaded prefix ++ and postfix ++ operators for             //
//   ProgCounter.                                                  //
// - added reset as a member function of ProgCounter.              //
// - added functionality to print reports to a file.               //
// 11/29/2017 - Brennan Young                                      //
// - added library dependencies: ctime, iostream.                  //
// 03/12/2019 - Brennan Young                                      //
// - constructor and reset now sets last time to 0, rather than    //
//   the time of the reset.                                        //
// 11/27/2019 - Brennan Young                                      //
// - added ProgCounter nesting with push() and pop().              //
// 12/02/2019 - Brennan Young                                      //
// - report returns address of the counter.                        //
// - log (file name) is now a public member.                       //
// - beginLog and continueLog no longer take a file name argument; //
//   they use the public log member.                               //
// - divided SUB into SUBP and SUBF.                               //
// - added SUBP to default mode.                                   //
// 12/03/2019 - Brennan Young                                      //
// - SUBP and SUBF now only aggregate progress up to levels where  //
//   the respective SUBP or SUBF flag is set.                      //
// - added elapsedTime method for determine which start time       //
//   should be used.                                               //
// 04/07/2020 - Brennan Young                                      //
// - corrected compiler warnings.                                  //
// 02/19/2021 - Brennan Young                                      //
// - push() can now be called to specify how much to change the    //
//   reporting level for the new nested ProgCounter object.        //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STDLIB_PROGRESSCOUNTER_20191126
#define YOUNG_STDLIB_PROGRESSCOUNTER_20191126

#include <ctime>    // time
#include <iostream> // std::cout
#include <fstream>  // std::fstream
#include <sstream>  // std::stringstream
#include <string>   // std:: string

namespace bystd { // Brennan Young standard namespace

// PROGRESS COUNTER OBJECt //////////////////////////////////////////

class ProgCounter {
private:
    // higher level (master counter)
    ProgCounter * master;
    bool masterFlag;
    
    // reporting time
    time_t last;
    
    // log file
    std::ofstream logf;
    
    // report string
    void leadingWhitespace (std::stringstream&) const;
    void ratioString (std::stringstream&) const;
    void timeString (std::stringstream&, const time_t&) const;
    time_t timeRemaining (const time_t&) const;
public:
    // bit flags for display mode (what to display in report)
    static const unsigned char PERC = 1;   // show percent flag
    static const unsigned char FRAC = 2;   // show fraction flag
    static const unsigned char TIME = 4;   // show time
    static const unsigned char TMIN = 8;   // T-minus (remaining)
    static const unsigned char SUBP = 16;  // hierarichal percent
    static const unsigned char SUBF = 32;  // hierarichal fraction
    static const unsigned char NL   = 64;  // new line
    static const unsigned char PLOG = 128; // print non-final to log
    
    // reporting
    unsigned char level;                // indentation level
    unsigned char reportLevel;          // maximum level to report
    
    // details
    std::string title;
    int begin, current, end;            // measurement for progress
    time_t start, interval;             // in seconds
    unsigned char mode;                 // bit flags
    std::string log;                    // log filename
    
    // constructors/destructor
    ProgCounter();                      // default constructor
    ProgCounter(const ProgCounter&);    // copy constructor
    ~ProgCounter ();                    // destructor
    
    // operators
    ProgCounter& operator= (const ProgCounter&);
    ProgCounter& operator++ ();
    ProgCounter operator++ (int);
    
    // management
    void reset();                       // reset to 0
    void pop ();                        // remove level
    void push ();                       // add new level
    void push (int);
    void beginLog ();                   // start empty log file
    void continueLog ();                // append to existing
    void endLog ();                     // close log file
    
    // report
    time_t elapsedTime () const;        // return elapsed time
    float fraction () const;            // return fraction complete
    ProgCounter& report ();
    void finalReport ();
}; // ProgCounter

// CONSTRUCTORS

// Default constructor
ProgCounter::ProgCounter ()
: master(NULL), masterFlag(false), last(time(NULL)),
    level(0), reportLevel(-1), title(""),
    begin(0), current(0), end(1),
    start(time(NULL)), interval(1),
    mode(ProgCounter::PERC | ProgCounter::SUBP | ProgCounter::TIME)
{}

// Copy constructor
ProgCounter::ProgCounter ( const ProgCounter & p )
: last(p.last), level(p.level), reportLevel(p.reportLevel),
    title(p.title), begin(p.begin), current(p.current), end(p.end),
    start(p.start), interval(p.interval), mode(p.mode), log(p.log)
{
    masterFlag = p.masterFlag;
    if ( masterFlag ) master = new ProgCounter(*p.master);
    else master = NULL;
    
    if ( p.logf.is_open() ) continueLog();
}

// Destructor
ProgCounter::~ProgCounter ()
{
    if ( masterFlag ) delete master;
    if ( logf.is_open() ) logf.close();
}

// OPERATORS

// Assignment operator.
ProgCounter & ProgCounter::operator= ( const ProgCounter & p )
{masterFlag = p.masterFlag;
    if ( masterFlag ) master = new ProgCounter(*p.master);
    else master = NULL;
    
    level       = p.level;
    reportLevel = p.reportLevel;
    begin       = p.begin;
    current     = p.current;
    end         = p.end;
    start       = p.start;
    last        = p.last;
    interval    = p.interval;
    title       = p.title;
    mode        = p.mode;
    log         = p.log;
    if ( p.logf.is_open() ) continueLog();
    
    return *this;
}

// prefix ++ (pre-increment)
ProgCounter & ProgCounter::operator++ ()
{ 
    ++current;
    return *this;
}

// postfix ++ (post-increment)
ProgCounter ProgCounter::operator++ ( int x )
{
    ProgCounter p = *this;
    ++current;
    return p;
}

// MANAGEMENT

void ProgCounter::reset()
{
    current = 0;
    start = time(NULL);
    last = 0;
}

void ProgCounter::pop ()
{
    if ( masterFlag ) {
        ProgCounter p = *master;
        p.last = last;
        delete master;
        *this = p;
    }
}

void ProgCounter::push ()
{
    master = new ProgCounter (*this);
    masterFlag = true;
    ++level;
    reset();
    last = master->last;
}

void ProgCounter::push ( int i )
{
    push();
    --level;
}

void ProgCounter::beginLog ()
{
    if ( logf.is_open() ) logf.close();
    if ( log.size() == 0 ) return;
    logf.open(log.c_str());
}

void ProgCounter::continueLog ()
{
    if ( logf.is_open() ) return;
    if ( log.size() == 0 ) return;
    logf.open(log.c_str(), std::ofstream::app);
}

void ProgCounter::endLog ()
{
    if ( logf.is_open() ) logf.close();
}

// REPORT STRING

// leading whitespace for report
void ProgCounter::leadingWhitespace ( std::stringstream & ss ) const
{
    for ( int i = 0; i < level; ++i ) ss << "  ";
}

// compute progress ratio string
void ProgCounter::ratioString ( std::stringstream & ss ) const
{
    const ProgCounter * p = this;
    bool flag = true;
    std::string temp;
    
    while ( flag ) {
        std::stringstream sstmp;
        sstmp << p->current << "/" << p->end;
        if ( temp.size() > 0 ) sstmp << "+" << temp;
        temp = sstmp.str();
        
        flag = p->masterFlag && (mode & SUBF);
        p = p->master;
    }
    
    ss << temp;
}

// determine time string
void ProgCounter::timeString ( std::stringstream & ss,
    const time_t & elapsed ) const
{
    if ( elapsed < 60 ) ss << elapsed << " seconds";
    else if ( elapsed < 3600 ) ss << (1.0 * elapsed / 60.0) << " minutes";
    else ss << (1.0 * elapsed / 3600.0) << " hours";
}

// estimate remaining time
time_t ProgCounter::timeRemaining ( const time_t & elapsed ) const
{
    float rate = 1.0 * elapsed / (current - begin);
    return (time_t) (rate * (end - current));
}

// REPORT

// compute elapsed time
time_t ProgCounter::elapsedTime () const
{
    const ProgCounter * p = this;
    time_t t0 = start;                          // starting time
    bool flag = p->masterFlag && (p->mode & (SUBF | SUBP));
    p = p->master;
    
    while ( flag ) {
        t0 = p->start;
        
        flag = p->masterFlag && (p->mode & (SUBF | SUBP));
        p = p->master;
    }
    
    return difftime(time(NULL), t0);
}

// compute fraction complete
float ProgCounter::fraction () const
{
    const ProgCounter * p = this;
    float f = 0.0f;                             // fraction complete
    bool flag = true;
    
    while ( flag ) {
        if ( p->end - p->begin > 0 ) {
            float x = p->end - p->begin;        // range
            float mf = (p->current - p->begin) / x;
            f = mf + f / x;                     // overall fraction
        }
        
        flag = p->masterFlag && (p->mode & SUBP);
        p = p->master;                          // go to master level
    }
    
    return f;
}

// print a progress report to the command window
ProgCounter & ProgCounter::report ()
{
    if ( level > reportLevel ) return *this;
    if ( difftime(time(NULL), last) < interval ) return *this;
    
    last = time(NULL);
    time_t elapsed = elapsedTime();
    int i;
    std::stringstream ss;
    
    ss << "\r";             // beginning of line
    leadingWhitespace(ss);  // indentation
    ss << title;            // title
    
    // time report, according to report mode
    if ( mode & FRAC ) {
        if ( title.size() > 0 ) ss << " ";
        ratioString(ss);
    }
    if ( mode & PERC ) ss << " " << (100.0f * fraction()) << "%";
    if ( mode & TIME ) {
        ss << " (";
        timeString(ss, elapsed);
        if ( mode & TMIN ) ss << "; ";
        else ss << ")";
    }
    if ( (mode & TMIN) && 0 < elapsed ) {
        if ( !(mode & TIME) ) ss << " (";
        ss << "est. ";
        timeString(ss, timeRemaining(elapsed));
        ss << " remaining)";
    }
    if ( (mode & PLOG) && logf.is_open() )
        logf << ss.str() << std::endl;
    
    // tailing whitespace
    int n = (80 - ss.str().length() );
    for ( i = 0; i < n; ++i ) ss << " ";
    for ( i = 0; i < n; ++i ) ss << "\b";
    
    if ( mode & NL ) ss << std::endl;
    
    std::cout << ss.str();
    
    return *this;
}

void ProgCounter::finalReport ()
{
    if ( level > reportLevel ) return;
    unsigned char pm = mode;
    time_t pi = interval;
    mode = mode | TIME | NL | PLOG;
    mode = mode & ~TMIN;
    interval = 0;
    report();
    mode = pm;
    interval = pi;
}

} // namespace bystd

#endif // YOUNG_STDLIB_PROGRESSCOUNTER_20191126