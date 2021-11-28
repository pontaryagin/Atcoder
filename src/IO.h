#include "MyHeader.h"

istream* pSTDIN = &cin;
ostream* pSTDOUT = &cout;

template<class T>
struct Reader
{
    static T _read(){
        T res;
        (*pSTDIN) >> res;
        return res;
    }
};

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

template<class T, class... Input>
T read(Input&&... args){
    return Reader<T>::_read(forward<Input>(args)...);
}

template<class T>
struct Writer
{
    static void _write(const T& arg){
        (*pSTDOUT) << arg;
    }
};

template<class T>
struct Writer<vector<T>>
{
    template<class... Del>
    static void _write(const vector<T>& arg, const string& del = " ", const Del&... dels){
        rep(i, 0, arg.size()){
            if (i != 0) (*pSTDOUT) << del;
            Writer<T>::_write(arg[i], dels...);
        }
    }
};

template<class T, class U>
struct Writer<pair<T, U>>
{
    template<class... Del1, class... Del2>
    static void _write(const pair<T, U>& arg, const string& del = "\n", tuple<Del1...> del1= tuple<>{}, tuple<Del2...> del2 = tuple<>{}){
        apply([&](auto &&... dels) { Writer<T>::_write(arg.first, dels...); }, del1);
        (*pSTDOUT) << del;
        apply([&](auto &&... dels) { Writer<U>::_write(arg.second, dels...); }, del2);
    }
    static void _write(const pair<T, U>& arg, const string& del, const string& del1, const string& del2){
        _write(arg, del, forward_as_tuple(del1), forward_as_tuple(del2));
    }
};

template<class T, class... Del>
void write(const T& arg, const Del&... del){
    Writer<T>::_write(arg, del...);
    (*pSTDOUT) << "\n";
}
