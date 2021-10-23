#include <type_traits>
#include <typeinfo>
#include <string>

template <typename T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::string r = typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}
//
//test(0);
//int x = 0;
//test(x);
//int&& y = 1;
//test(y);
//const int a = 0;
//test(a);
//const int& b = 1;
//test(b);
//const int&& v = 1;
//test(v);
