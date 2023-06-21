#include <iostream>
#include <time.h>
#include <cstdlib>
#include <locale.h>
#include <chrono>
using namespace std;
typedef unsigned char BASE;
typedef unsigned short int DBASE;
typedef unsigned int QBASE;
#define BASE_SIZE (sizeof(BASE)*8);
#define DBASE_SIZE (sizeof(DBASE)*8);
#define QBASE_SIZE (sizeof(DBASE)*8);

class BN
{
	BASE* coef;
	int len, maxlen;
public:
	BN(int nmax, int t);
	BN(const BN&);
	~BN() { delete[] coef; };
	BN operator =(const BN&);
	bool operator ==(const BN&);
	bool operator !=(const BN&);
	bool operator >(const BN&);
	bool operator >=(const BN&);
	bool operator <(const BN&);
	bool operator <=(const BN&);
	BN operator +(const BASE&);
	BN operator +(const BN&);
	BN operator -(const BN&);
	BN operator *(const BASE&);
	BN operator *(BN&);
	BN operator /(const BASE&);
	BN operator /(BN);
	BN operator %(BN);
	BASE operator %(const BASE&);
	BN& operator +=(const BN&);
	BN& operator -=(const BN&);
	BN Squaring();
	BN Pow(int);
	BN Mod(BN);
	int GetLen() { return len; };
	BN Obr1(int);
	BN Obr2(int);
	void Output_16();
	void Input_16();
	bool MillerRabin(int);
	BN ModExp(BN, BN);
	BN phifun();
	bool Ferma(int);
	int NumberOfDigit(BASE);
	friend int test();
	BN Exp(BN);
	BN Gordon();
	friend ostream& operator << (ostream&, BN&);
	friend istream& operator >> (istream&, BN&);

};

BN::BN(int nmax = 1, int t = 1)
{
	maxlen = nmax;
	len = maxlen;
	int i, k = 0, n = BASE_SIZE;
	coef = new BASE[maxlen];
	if (t == 1) for (i = 0; i < maxlen; i++) coef[i] = 0;
	if (t == 2)
	{
		for (i = 0; i < maxlen; i++)
		{
			coef[i] = 0;
			k = 0;
			while (k < n)
			{
				coef[i] |= rand() << k;
				k += 16;
			}
		}
	}
	for (i = maxlen - 1; i > 0 && coef[i] == 0; i--) len--;
}

BN::BN(const BN& bn)
{
	maxlen = bn.maxlen;
	len = bn.len;
	coef = new BASE[maxlen];
	for (int i = len - 1; i >= 0; i--) coef[i] = bn.coef[i];
}

BN BN::operator=(const BN& bn)
{
	if (this != &bn) {
		delete[] coef;
		maxlen = bn.maxlen;
		len = bn.len;
		coef = new BASE[maxlen];
		for (int i = len - 1; i >= 0; i--) coef[i] = bn.coef[i];
	}
	return *this;
}

bool BN::operator==(const BN& b)
{
	if (len != b.len) return false;
	for (int i = len - 1; i >= 0; i--) if (coef[i] != b.coef[i]) return false;
	return true;
}

bool BN::operator!=(const BN& b)
{
	if (len != b.len) return true;
	for (int i = len - 1; i >= 0; i--) if (coef[i] != b.coef[i]) return true;
	return false;
}

bool BN::operator>=(const BN& b)
{
	if (len < b.len) return false;
	if (len > b.len) return true;
	for (int i = len - 1; i >= 0; i--) {
		if (coef[i] > b.coef[i]) return true;
		if (coef[i] < b.coef[i]) return false;
	}
	return true;
}

bool BN::operator>(const BN& b)
{
	if (len < b.len) return false;
	if (len > b.len) return true;
	for (int i = len - 1; i >= 0; i--) {
		if (coef[i] > b.coef[i]) return true;
		if (coef[i] < b.coef[i]) return false;
	}
	return false;
}

bool BN::operator<=(const BN& b)
{
	if (len > b.len) return false;
	if (len < b.len) return true;
	for (int i = len - 1; i >= 0; i--) {
		if (coef[i] < b.coef[i]) return true;
		if (coef[i] > b.coef[i]) return false;
	}
	return true;
}



bool BN::operator<(const BN& b)
{
	if (len > b.len) return false;
	if (len < b.len) return true;
	for (int i = len - 1; i >= 0; i--) {
		if (coef[i] < b.coef[i]) return true;
		if (coef[i] > b.coef[i]) return false;
	}
	return false;
}

BN BN::operator+(const BASE& v)
{
	int j = 0, k = 0, n = BASE_SIZE;
	BN w(len + 1, 1);
	w.len = w.maxlen;
	DBASE tmp;
	tmp = (DBASE)coef[j] + v;
	w.coef[j] = (BASE)tmp;
	k = (BASE)(tmp >> n);
	j++;
	while (j < len) {
		tmp = (DBASE)coef[j] + k;
		w.coef[j] = (BASE)tmp;
		k = (BASE)(tmp >> n);
		j++;
	}
	w.coef[j] = k;
	for (j = w.maxlen - 1; j > 0 && w.coef[j] == 0; j--) w.len--;
	return w;
}

BN BN::operator+(const BN& v)
{
	int k = 0, j = 0, l = max(len, v.len) + 1, t = min(len, v.len), n = BASE_SIZE;
	BN w(l, 1);
	w.len = w.maxlen;
	DBASE tmp;
	while (j < t)
	{
		tmp = (DBASE)coef[j] + (DBASE)v.coef[j] + k;
		w.coef[j] = (BASE)tmp;
		k = (BASE)(tmp >> n);
		j++;
	}
	while (j < len)
	{
		tmp = (DBASE)coef[j] + k;
		w.coef[j] = (BASE)tmp;
		k = (BASE)(tmp >> n);
		j++;
	}
	while (j <
		v.len)
	{
		tmp = (DBASE)v.coef[j] + k;
		w.coef[j] = (BASE)tmp;
		k = (BASE)(tmp >> n);
		j++;
	}
	w.coef[j] = k;
	if (k == 0) w.len--;
	return w;
}

BN BN::operator-(const BN& v)
{
	if (*this < v)
	{
		cout << "Ошибка - второе число больше превого" << endl;
		cout << *this << "-" <<endl;
		exit(0);
	}
	int j = 0, k = 0, n = BASE_SIZE;
	DBASE tmp;
	BN w(len, 1);
	w.len = w.maxlen;
	while (j < v.len)
	{
		tmp = ((DBASE)1 << n) | coef[j];
		tmp = tmp - (DBASE)v.coef[j] - k;
		w.coef[j] = (BASE)tmp;
		k = !(tmp >> n);
		j++;
	}
	while (j < len)
	{
		tmp = ((DBASE)1 << n) | coef[j];
		tmp -= k;
		w.coef[j] = (BASE)tmp;
		k = !(tmp >> n);
		j++;
	}
	for (j = w.maxlen - 1; j > 0 && w.coef[j] == 0; j--) w.len--;
	return w;
}

BN BN::operator*(const BASE& v)
{
	BN w(len + 1, 1);
	w.len = w.maxlen;
	int j, n = BASE_SIZE; BASE k = 0;
	DBASE tmp;
	for (j = 0; j < len; j++)
	{
		tmp = (DBASE)coef[j] * v + k;
		w.coef[j] = (BASE)tmp;
		k = (BASE)(tmp >> n);
	}
	w.coef[j] = k;
	if (k == 0) w.len--;
	return w;
}

BN BN::operator*(BN& v)
{
	BN w(len + v.len, 1), z;
	w.len = w.maxlen;
	int j, i, n = BASE_SIZE;
	BASE k = 0;
	DBASE tmp;
	if (*this == z || v == z) {
		w = z; return w;
	}
	for (j = 0; j < v.len; j++) {
		k = 0;

		for (i = 0; i < len; i++) {
			tmp = (DBASE)coef[i] * (DBASE)v.coef[j] + (DBASE)w.coef[i + j] + k;
			w.coef[i + j] = (BASE)tmp;
			k = (BASE)(tmp >> n);
		}
		w.coef[j + len] = k;
	}
	if (k == 0) w.len--;
	return w;
}

BN BN::operator/(const BASE& v)
{
	BN q(len, 1);
	q.len = q.maxlen;
	int j, n = BASE_SIZE;
	BASE r = 0;
	DBASE tmp;
	for (j = 0; j < len; j++) {
		tmp = ((DBASE)r << n) | (DBASE)coef[len - 1 - j];
		q.coef[len - 1 - j] = (BASE)(tmp / v);
		r = (BASE)(tmp % v);
	}
	for (j = q.maxlen - 1; j > 0 && q.coef[j] == 0; j--) q.len--;
	return q;
}

BN BN::operator/(BN v)

{
	BN r(v.len, 1), w;
	if (*this < v) {
		r = *this; return w;
	}
	if (*this == v) {
		r = w; w = w + 1; return w;
	}
	if (v == w) {
		cout << "Ошибка - деление на ноль"; exit(0);
	}
	int m = len - v.len, n = BASE_SIZE;
	BASE d; DBASE b = pow(2, n);
	BN q(m + 1, 1), kop(*this);

	if (v.len == 1) {
		BASE del = v.coef[0];
		q = kop / del;
		return q;
	}

	q.len = q.maxlen; r.len = r.maxlen;
	d = (b / ((double)v.coef[v.len - 1] + 1));
	int num = kop.len;
	kop = kop * d; v = v * d;
	if (kop.len == num) {
		kop.len++; kop.coef[kop.len - 1] = 0;
	}

	for (int j = m; j >= 0; j--) {
		DBASE q1; DBASE r1;
		q1 = ((DBASE)kop.coef[j + v.len] * b + (DBASE)kop.coef[j + v.len - 1]) / (DBASE)v.coef[v.len - 1];
		r1 = ((DBASE)kop.coef[j + v.len] * b + (DBASE)kop.coef[j + v.len - 1]) % (DBASE)v.coef[v.len - 1];
		if (q1 == b || q1 * v.coef[v.len - 2] > b * r1 + kop.coef[j + v.len - 2]) {
			q1--; r1 = r1 + v.coef[v.len - 1];
		}
		if (r1 < b && (q1 == b || q1 * v.coef[v.len - 2] > b * r1 + kop.coef[j + v.len - 2])) {
			q1--; r1 = r1 + v.coef[v.len - 1];
		}
		BN u1(v.len + 1, 1);
		u1.len = u1.maxlen;
		int k = v.len;
		int i;
		for (i = 0; i <= k; i++) {
			u1.coef[i] = kop.coef[j + i];
		}
		if (u1 < v * q1) { q1--; }
		u1 -= v * q1;
		for (i = 0; i <= k; i++) {
			kop.coef[j + i] = u1.coef[i];
		}

		for (i = u1.len; i <= v.len; i++) {
			kop.coef[j + i] = 0;
		}
		q.coef[j] = q1;
	}

	return q;
}

BN BN::operator%(BN v)
{
	int m = len - v.len, n = BASE_SIZE;
	BASE d; DBASE b = pow(2, n);
	BN r(v.len, 1), kop(*this), z;

	if (v == z) {
		cout << "Ошибка - деление на ноль"; exit(0);
	}

	if (*this < v) {
		r = *this; return r;
	}
	if (*this == v) {

		r = z; return r;
	}
	if (v.len == 1) {
		BASE del = v.coef[0];
		del = kop % del;
		r.coef[0] = del;
		r.len = 1;
		return r;
	}
	r.len = r.maxlen;
	d = (b / ((double)v.coef[v.len - 1] + 1));
	int num = kop.len;
	kop = kop * d; v = v * d;
	if (kop.len == num) {
		kop.len++; kop.coef[kop.len - 1] = 0;
	}
	for (int j = m; j >= 0; j--) {
		DBASE q1; DBASE r1;
		q1 = ((DBASE)kop.coef[j + v.len] * b + (DBASE)kop.coef[j +
			v.len - 1]) / (DBASE)v.coef[v.len - 1];
		r1 = ((DBASE)kop.coef[j + v.len] * b + (DBASE)kop.coef[j + v.len - 1]) % (DBASE)v.coef[v.len - 1];
		if (q1 == b || q1 * v.coef[v.len - 2] > b * r1 + kop.coef[j + v.len - 2]) {
			q1--; r1 = r1 + v.coef[v.len - 1];
		}
		if (r1 < b && (q1 == b || q1 * v.coef[v.len - 2] > b * r1 + kop.coef[j + v.len - 2])) {
			q1--; r1 = r1 + v.coef[v.len - 1];
		}
		BN u1(v.len + 1, 1);
		u1.len = u1.maxlen;
		int k = v.len;
		int i;
		for (i = 0; i <= k; i++) {
			u1.coef[i] = kop.coef[j + i];
		}
		if (u1 < v * q1) { q1--; }
		u1 -= v * q1;
		for (i = 0; i <= k; i++) {
			kop.coef[j + i] = u1.coef[i];
		}

		for (i = u1.len; i <= v.len; i++) {
			kop.coef[j + i] = 0;
		}
	}
	kop = kop / d;
	return kop;
}

BASE BN::operator%(const BASE& v)
{
	BN q(len, 1);
	q.len = q.maxlen;
	int j, n = BASE_SIZE;
	BASE r = 0;
	DBASE tmp;
	for (j = 0; j < len; j++) {
		tmp = ((DBASE)r << n) | (DBASE)coef[len - 1 - j];
		r = (BASE)(tmp % v);
	}
	return r;
}


BN& BN::operator+=(const BN& b)
{
	*this = *this + b;
	return *this;
}

BN& BN::operator-=(const BN& b)
{
	*this = *this - b;
	return *this;
}

BN BN::Squaring()
{
	BN y(len*2+1, 1);
	y.len = y.maxlen;
	for (int i = 0; i <= 2 * len - 1; i++)
	{
		y.coef[i] = 0;
	}
	
	DBASE uv;
	QBASE cuv;
	int shift = BASE_SIZE; 
	
	for (int i = 0; i <= len - 1; i++)
	{
		uv = (DBASE)y.coef[2 * i] + (DBASE)coef[i] * (DBASE)coef[i];
		y.coef[2 * i] = (BASE)uv;
		cuv =(QBASE)uv;
		for (int j = i + 1; j <= len - 1; j++)
		{
			cuv = (QBASE)y.coef[i + j] + 2 * (QBASE)coef[i] * (QBASE)coef[j] + (QBASE)(cuv >> shift);
			y.coef[i + j] = (BASE)(cuv);

		}
		uv = (DBASE)(cuv >> shift);
		y.coef[i + len + 1] += (BASE)(uv >> shift);
		
		y.coef[i + len] += (BASE)uv;
	}
	int i = y.len - 1;
	while (y.coef[i] == 0)
	{
		y.len--; i--;
	}
	return y;
}

BN BN::Pow(int y)
{
	
	BN z = *this;

	int j = 1;
	int i = y;
	while (i > 1)
	{
		i=i>>1;
		j++;
	}
	//cout << j<<endl;
	
	for (i = j - 2; i >= 0; i--)
	{
		z = z.Squaring();
		if ((((1<<i)&y)>>i)== 1)z = z * (*this);
	}
	return z;
}

BN BN::Mod(BN m)
{
	BN x = *this;

	//cout << "x: " << x << endl;
	
	// во второй коэффициент занести 1
	if ((m.coef[m.len - 1]==0)||(len>2*m.len))
	{
		exit(0);
	}

	BN z(2, 1);
	z.len = z.maxlen;
	BASE ones = 1;
	z.coef[1] = ones;//для получения основания системы счисления
	z = z.Pow(2*m.len)/m;
	
	BN q1;
	q1 = (x.Obr1(m.len - 1)*z).Obr1(m.len+1);
	//cout <<q1 << endl;
	BN r1(len,1);
	BN r2(len,1);
	//cout << "x: "<<x<<endl;
	r1 = x.Obr2(m.len + 1);
	r2 = (q1 * m).Obr2(m.len + 1);
	//cout << r1 << endl;
	//cout << r2 << endl;
	BN r(2,1);
	r.len = r.maxlen;
	if (r1 >= r2)
	{
		r = r1 - r2;
	}
	else
	{
		r.coef[1] = ones;
		//cout << 'a';
		r = r.Pow(m.len+1)+r1-r2;
	}
	//cout <<endl<< r<<endl;

	while (r >= m)
	{
		r = r - m;
	}

	return r;

}

BN BN::Obr1(int i)// деление на b^i
{
	BN c(1,1);
	if (i >= len)return c;
	BN b(len - i, 1);
	b.len = b.maxlen;
	for (int j = len - 1; j >= i; j--) { b.coef[j - i] = coef[j];};
	//b.Output_16();
	return b;
}

BN BN::Obr2(int i)//mod b^i
{
	BASE b = 0;
	BN a(i, 1);
	a.len = a.maxlen;
	for (int j = 0; j < i && j<len; j++) { a.coef[j] = coef[j]; };
	return a;
}

BN BN::ModExp(BN y, BN mod)
{
	BN z(*this);
	int base = BASE_SIZE;
	int n = base * (y.len - 1), n1 = -2;
	BASE y_coef = y.coef[y.len - 1];

	while (y_coef != 0)
	{
		y_coef >>= 1;
		n++; n1++;
	}
	y_coef = y.coef[y.len - 1];
	for (int i = n - 2, j = y.len - 1; i >= 0; i--)
	{
		if (n1 < 0)
		{
			n1 = base - 1;
			y_coef = y.coef[--j];

		}
		z = z * z; z = z.Mod(mod);

		
		if (((1 << n1) & y_coef) != 0)
		{
			z = z * (*this); z = z.Mod(mod);
			
		}
		n1--;
	}
	if (z > mod) z = z.Mod(mod);
	return(z);
}

int BN :: NumberOfDigit(BASE y)
{
	int n = 0;
	while (y != 0)
	{
		y = y / 2;
		n++;
	}
	return n;
}

BN BN::Exp(BN y)
{
	BN z; int base = BASE_SIZE;
	int n = base * (y.len - 1) + NumberOfDigit(y.coef[y.len - 1]), mask = 0;
	if (NumberOfDigit(y.coef[y.len - 1]) - 2 > 0) mask = 1 << NumberOfDigit(y.coef[y.len - 1]) - 2;
	z = *this;
	for (int i = n - 2, j = y.len - 1; i >= 0; i--)
	{
		if (mask == 0)
		{
			if (j > 0) j--; mask = 1 << base - 1;
		}
		z = z.Squaring();
		int b = mask & y.coef[j];
		if (mask & y.coef[j]) z = z * (*this);
		mask >>= 1;

	}
	return z;
}

BN BN::phifun()
{
	BN result = *this, n = *this, z, i, j, b;
	i.coef[0] = 2;
	j.coef[0] = 1;
	z.coef[0] = 0;
	for (i; i.Squaring() <= n; i += j)
	{
		BN t = i.Squaring();
		//cout << t << endl;
		if (n % i == z)
		{
			while (n % i == z)
				n = n / i;
			result = result - result / i;
		}
	}

	if (n > j)
	{
		result = result - result / n;
	}
	return result;
}

bool BN::Ferma(int t)
{
	int kof = t;
	BN b(len, 1);
	BN r(len, 1);
	b.coef[0] = 0; 
	r.coef[0] = 2;
	if (((*this) % r == b) && (*this != r))
	{ 
		cout << "Число составное" << endl;
		return false;
	}
	b.coef[0] = 3;
	if (*this <= b) 
	{
		return false;
	}
	for (t; t > 0; t--)
	{
		BN a(len, 2); 
		b.coef[0] = 2;
		r.coef[0] = 3;
		a = b + a.Mod(*this - r);//2<=a<=n-2
		b.coef[0] = 1;
		r = a.ModExp(*this - b, *this);//r=a^(n-1)mod n
		if (r != b) 
		{
			cout << "Число составное" << endl;
			return false;
		}
	}
	//BN phi = phifun().Pow(kof), n = Pow(kof);
	cout << "Число, вероятно, простое" << endl;
	//cout << "Вероятность ошибки: " << phi << "/" << n << endl;
	return true;
}

bool BN::MillerRabin(int t)
{
	int kof = t;
	int s = 1; BN r, с, d, e;
	с.coef[0] = 1;
	d.coef[0] = 2;
	e.coef[0] = 0;
	if (*this <= d + с) //n<=3
	{ 
		//cout<<"Алгоритм не применим для заданного числа";
		return false;
	}

	if (*this % d == e) //если чётный
	{ 
		//cout << "Составное";
		return false; 
	}

	r = *this / d;
	while (r % d != с)//n=2^s*r+1
	{
		r = r / d;
		s++;
	}

	for (int i = 0; i < t; i++)
	{
		BN b(1, len);
		b = b.Mod(*this - d - d) + d;//2<=b<=n-2
		BN y = b.ModExp(r, *this);
		if (y != с && y != *this - с)//y!=1, y!=-1 mod n
		{
			int j = 1;
			while (j < s && y != *this - с)
			{
				y = y.ModExp(d, *this);
				if (y == с) 
				{ 
					//cout << "Составное";
					return false; 
				}
				j++;
			}
			if (y != *this - с) 
			{ 
				//cout << "Составное";
				return false;
			}
		}
	}
	/*BN a;
	a.coef[0] = 4;
	BN phi = phifun().Pow(kof), n = Pow(kof) * a;*/
	//cout << "Число, вероятно, простое" << endl;
	//cout << "Вероятность ошибки: " << phi<<'/'<<n<<endl;
	return true;
}

BN BN::Gordon()
{
	BN s(5, 2);
	BN t(5, 2);
	while (!s.MillerRabin(10))//генерация простого числа
	{
		BN a(5, 2);
		s = a;
	}
	while (!t.MillerRabin(10))//генерация ещё одного
	{
		BN a(5, 2);
		t = a;
	}
	BN i0(1, 2);//генерация целого случайного числа
	BN r;
	BN p0;
	BN f;
	BN h;
	f.coef[0] = 2;
	h.coef[0] = 1;
	while (!((f * i0 * t + 1).MillerRabin(10)))//ищем i(i0),при котором 2it+1 станет простым 
	{
		i0 += h;
	}
	r = f * i0 * t + 1;//считаем простое r
	p0 = s.ModExp(r - f, r) * f * s - h;
	BN j0(1, 2);
	while (!((f * j0 * r * s + p0).MillerRabin(10)))
	{
		j0 += h;
	}
	p0 = f * j0 * r * s + p0;
	f.coef[0] = 0;
	//cout << p0 << endl;
	if (((p0 - h) % r == f) && ((p0 + h) % s == f) && ((r - h) % t == f))
	{
		cout << "Cильное простое: " << endl;
	}
	return p0;
}

void BN::Output_16()
{
	int j, i = 0, tmp, n = BASE_SIZE;
	int k = BASE_SIZE; k -= 4;
	char* s = new char[len * (n / 4)];
	for (j = len - 1; j >= 0;) {
		tmp = (coef[j] >> k) & (0xf);
		if (tmp >= 0 && tmp <= 9) s[i] = (char)(tmp + '0');
		if (tmp >= 10 && tmp <= 15) s[i] = (char)(tmp - 10 + 'a');
		i++; k -= 4;
		if (k < 0) {
			k = BASE_SIZE; k -= 4; j--;
		}
	}
	for (i = 0; s[i] == '0' && i < n / 4 - 1; i++);
	for (i; i < len * (n / 4); i++) cout << s[i];
	cout << endl;
}

void BN::Input_16()
{
	int i, j = 0, k = 0, tmp = 0, n = BASE_SIZE;
	char buf[256];
	delete[] coef;
	cout << "Введите число в 16-ой системе" << endl;
	gets_s(buf);
	len = (strlen(buf) - 1) / (n / 4) + 1;
	maxlen = len;
	coef = new BASE[maxlen];
	for (i = 0; i < maxlen; i++) coef[i] = 0;
	for (i = strlen(buf) - 1; i >= 0; i--)
	{
		if (buf[i] >= '0' && buf[i] <= '9') tmp = buf[i] - '0';
		if (buf[i] >= 'a' && buf[i] <= 'f') tmp = buf[i] - 'a' + 10;
		if (buf[i] >= 'A' && buf[i] <= 'F') tmp = buf[i] - 'A' + 10;
		if (buf[i] < '0' || (buf[i] > '9' && buf[i] < 'A') || (buf[i] > 'F' && buf[i] < 'a') || buf[i]>'f')
		{
			cout << "Ошибка,неверный символ" << endl; exit(0);
		}
		coef[j] |= tmp <<(k * 4);
		k++;
		if (k >= n / 4)
		{
			k = 0; j++;
		}
	}
	for (i = maxlen - 1; i > 0 && coef[i] == 0; i--) len--;
}

int test()
{
	int M = 1000, T = 1000, m, n;
	BN A, B, C, D;
	do 
	{
		n = rand() % M + 1; m = rand() % M + 1;
		BN A1(n, 2), B1(m, 2);
		A = A1; B = B1;
		C = A / B; D = A % B;
		cout << T << ' ' << m << ' ' << n << endl;
	} 
	while (A == C * B + D && A - D == C * B && D < B && --T);
	return T;
}

istream& operator >>(istream& p, BN& u)
{
	BN v;
	int j; BASE c = 10, t;
	char buf[256];
	cout << "\nВведите число" << endl;
	gets_s(buf);
	for (j = 0; j < strlen(buf); j++)
	{
		t = buf[j] - '0';
		v = v * c + t;
	}
	u = v;
	return p;
}

ostream& operator <<(ostream& p, BN& u)
{
	BASE t, c = 10; int j = 0, i = 0; BN d, a(u);
	int n = BASE_SIZE;
	char* s = new char[(n + 1) / 3 * u.len];
	if (a == d) { cout << '0'; return p; }
	while (a != d) {
		t = a % c;
		s[j] = t + '0';
		a = a / c;
		j++;
	}
	for (i = j - 1; i >= 0; i--) cout << s[i];
	delete[] s;
	return p;
}


int main()
{
	setlocale(LC_ALL, "Russian");
	srand(time(NULL));
	
	/*BN blc;
	BN blb;
	cin >> blc;
	auto start = chrono::steady_clock::now();
	blb = blc.Squaring();
	auto end = chrono::steady_clock::now();
	auto result = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	cout << blb<<' '<<result.count()<<endl;
	start = chrono::steady_clock::now();
	blb = blc * blc;
	end = chrono::steady_clock::now();
	result = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	cout << blb << ' ' << result.count() << endl;*/
	
	/*BN a;
	BN b;
	int i;
	cin >> a;
	cin >> i;
	b = a;
	auto start = chrono::steady_clock::now();
	a=a.Pow(i);
	auto end = chrono::steady_clock::now();
	auto result = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	cout << a << ' ' << result.count() << endl;
	a = b;
	start = chrono::steady_clock::now();
	for (int j = 0; j < i-1; j++)
	{
		a = a * b;
	}
	end = chrono::steady_clock::now();
	result = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	cout << a << ' ' << result.count() << endl;*/

	/*BN a;
	BN b;
	BN c;
	BN a1;
	BN b1;
	BN c1;
	BN m;
	cin >> a;
	cin >> b;
	cin >> c;
	a1 = a;
	b1 = b;
	c1 = c;
	int k;
	cin >> m;
	BN base;
	cin >> base;
	base = base.Pow(2 * m.GetLen()) / m;

	auto start = chrono::steady_clock::now();
	a = a.Mod(m,base);
	b = b.Mod(m, base);
	c = c.Mod(m, base);
	auto end = chrono::steady_clock::now();
	auto result = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	cout << "Алгоритм Барретта a: " << a << endl;
	cout << "Алгоритм Барретта b: " << b << endl;
	cout << "Алгоритм Барретта c: " << c << endl;
	cout << endl;
	cout << "Время: " << result.count() << endl;
	cout << endl<<endl;
	
	start = chrono::steady_clock::now();
	a1 = a1 % m;
	b1 = b1 % m;
	c1 = c1 % m;
	end = chrono::steady_clock::now();
	result = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	cout << " a: " << a1 << endl;
	cout << " b: " << b1 << endl;
	cout << " c: " << c1 << endl;
	cout << endl;
	cout << "Время: " << result.count() << endl;*/

	BN k;
	BN n = k.Gordon();
	/*if (n.MillerRabin(10))
	{
		cout << "Простое"<<endl;
	}
	else
	{
		cout << "Составное";
	}*/
	cout << n;
}
