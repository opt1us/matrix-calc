#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<math.h>
#include<string>
#include<algorithm>
#include<numeric>


#define vi  vector<int>
#define vd  vector<double>
#define vpi vector<pair<int, int>>
#define pii pair<int, int>

using namespace std;

// Условие для всех задач - вывод должен полностью соответствовать стандартному решению в тетради
// или помогать его воспроизвести

// A1 Поиск всех подмножеств (подкортежей) из заданного множества с определенной суммой
// A2 Нахождение алфавитного дерева минимальной стоимости (задача расставления скобок, динамическое программирование)
// A3 Три жадных алгоритма для решения задачи о рюкзаке (по макс. цене, по мин. весу, по макс. удельной цене)
// A4 Решение задачи о рюкзаке при помощи динамического программирования

// Для задачи A1
vi vg(1, -1);

template<typename T>
void printv(const vector<T>& v, int pad) {
    for (T e: v) {
        cout << setw(pad) << e << ' ';
    }
    cout << '\n';
}

template<typename T>
void printv(const vector<T>& v) {
    for (T e: v) {
        cout << e << ' ';
    }
    cout << '\n';
}

template<typename T, typename P>
void printvp(const vector<pair<T, P>>& v) {
    for (auto e: v) {
        cout << e.first << ':' << e.second << ' ';
    }
    cout << '\n';
}

void recur(const vi& v, int i, int s, const int& d, vi vt, int vgindx) {
    vg[vgindx] = s;
    if (s == d) {
        printv(vt);
        cout << s << '\n';
        return;
    }
    if (i == v.size()) return;
    int sum = accumulate(v.begin() + i, v.end(), s);
    if (s + sum < d || s + v[i] > d) return;
    vt.push_back(i);
    recur(v, i + 1, s + v[i], d, vt, vgindx * 2);
    recur(v, i + 1, s, d, vi(vt.begin(), vt.end() - 1), vgindx * 2 + 1);
}

void solve_A1() {
    int n;
    int d;
    ifstream in("2.txt");
    in >> n;
    vi v(n);
    for (int i = 0; i < n; i++) {
        in >> v[i];
    }
    in >> d;
    printv(v);
    cout << d << '\n';
    sort(v.begin(), v.end());
    // solving
    ofstream out("bt.txt");
    cout << "Backtracking solves all subsets sum equal d\n\n";
    vi vt;
    vg.resize(1<<(n+1), -1);
    recur(v, 0, 0, d, vt, 1);
    for (int i = 0; i < n + 1; i++) {
        for (int j = (1<<i); j < (1<<(i+1)); j++) {
            int pad = (1<<(n + 2 - i));
            if (j == (1<<i)) pad >>= 1;
            if (vg[j] != -1) out << setw(pad) << vg[j];
            else out << setw(pad) << ' ';
        }
        out << '\n';
    }
}

void solve_A3() {
    int n;
    int d;
    ifstream in("2.txt");
    in >> n;
    vpi v(n);
    vpi w(n);
    for (int i = 0; i < n; i++) {
        in >> v[i].second >> w[i].second;
        v[i].first = w[i].first = i;
    }
    in >> d;
    printvp(v);
    printvp(w);
    cout << d << '\n';
    // solving max value
    {
        vpi v_(v);
        vpi w_(w);
        cout << "Greedy backpack (max value) #1\n\n";
        sort(v_.begin(), v_.end(), [&](pair<int, int> p1, pair<int, int> p2) {
            // if (p1.second == p2.second) return w_[p1.first].second < w_[p2.first].second;
            return p1.second > p2.second;
        });
        int dc = d;
        int s = 0;
        vi ei;
        for (int i = 0; i < n; i++) {
            int value = v_[i].second;
            int weight = w_[v_[i].first].second;
            if (weight > dc) continue;
            s += value;
            dc -= weight;
            ei.push_back(v_[i].first);
        }
        cout << "\tItems indices: "; printv(ei);
        cout << "\tFull Cost: " << s << '\n';
        cout << "\tTotal weight: " << d - dc << '\n';

    }
    cout << '\n';
    // min weight
    {
        vpi v_(v);
        vpi w_(w);
        cout << "Greedy backpack (min weight) #2\n\n";
        sort(w_.begin(), w_.end(), [&](pair<int, int> p1, pair<int, int> p2) {
            // if (p1.second == p2.second) return v_[p1.first].second > v_[p2.first].second;
            return p1.second < p2.second;
        });
        int dc = d;
        int s = 0;
        vi ei;
        for (int i = 0; i < n; i++) {
            if (w_[i].second > dc) continue;
            s += v_[w_[i].first].second;
            dc -= w_[i].second;
            ei.push_back(w_[i].first);
        }
        cout << "\tItems indices: "; printv(ei);
        cout << "\tFull Cost: " << s << '\n';
        cout << "\tTotal weight: " << d - dc << '\n';
    }
    cout << '\n';
    // max unit value
    {
        cout << "Greedy backpack (max unit weight) #3\n\n";
        vpi uv(n);
        for (int i = 0; i < n; i++) {
            uv[i].first = i;
            uv[i].second = v[i].second / w[i].second;
        }
        sort(uv.begin(), uv.end(), [](pair<int, int> p1, pair<int, int> p2) {
            return p1.second > p2.second;
        });
        int dc = d;
        int s = 0;
        vi ei;
        for (int i = 0; i < n; i++) {
            int idx = uv[i].first;
            if (w[idx].second > dc) continue;
            s += v[idx].second;
            dc -= w[idx].second;
            ei.push_back(idx);
        }
        cout << "\tItems indices: "; printv(ei);
        cout << "\tFull Cost: " << s << '\n';
        cout << "\tTotal weight: " << d - dc << '\n';
    }
}

int main() {
    solve_A1();
}
