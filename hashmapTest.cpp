#include <unordered_set>
#include <iostream>
using namespace std;

int main()
{
    unordered_set<int> S;
    int n, m, t;
    cin >> n;
    for (int i = 0; i < n; i++)
    {
        cin >> t;
        S.insert(t);
    }
    cin >> m;
    for (int i = 0; i < m; i++)
    {
        cin >> t;
        if (S.find(t) != S.end())
            cout << "YES" << endl;
        else
            cout << "NO" << endl;
    }
    cout << 1;
    return 0;
}