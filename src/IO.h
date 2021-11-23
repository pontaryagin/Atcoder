#include "MyHeader.h"

istream* pSTDIN = &cin;

template<class T, class Enabler = void>
struct Reader{};

template<class T>
struct Reader<vector<T>>
{
    template<class ...Input>
    static vector<T> _read(ll n, Input&&... args){
        vector<T> res(n);
        rep(i,0,n){
            res[i] = move(Reader<T>::_read(std::forward<Input>(args)...));
        }
        return res;
    }
};

template<class T, class U>
struct Reader<pair<T, U>>
{
    template<class... Input1, class... Input2>
    static pair<T, U> _read(tuple<Input1...> args1 = tuple<>{}, tuple<Input2...> args2 = tuple<>{}){
        auto l = apply([](auto &&... args) { return Reader<T>::_read(move(args)...); }, args1);
        auto r = apply([](auto &&... args) { return Reader<U>::_read(move(args)...); }, args2);
        return make_pair(move(l), move(r));
    }
    template<class Input1, class Input2>
    static pair<T, U> _read(Input1&& arg1, Input2&& arg2){
        return _read(forward_as_tuple(arg1), forward_as_tuple(arg2));
    }

};

template<class T>
struct Reader<T>
{
    static T _read(){
        T res;
        (*pSTDIN) >> res;
        return res;
    }
};

template<class T, class... Input>
T read(Input&&... args){
    return Reader<T>::_read(forward<Input>(args)...);
}




