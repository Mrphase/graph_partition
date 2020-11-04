// #pragma once
// #pragma once
// #pragma once
#include<iostream>
using namespace std;
template <typename Type>//模板函数，有效区域该声明下面一行代码
class Stack {
private:
	Type* stk; //起始地址
	int top; //始终指向栈顶元素
	int MAXN;//栈的最大存储容量
public:
	Stack(int size); //构造函数，初始化一个栈时需要指定初始大小
	~Stack();//析构函数
	int push(Type x);
	int pop(Type* px);
	Type getTop();
	int isEmpty()const;
	int isFull()const;
	int size()const;
};
template <typename Type>
//注意要加上尖括号模板，和普通变量类型不同
Stack<Type>::Stack(int size) {
	MAXN = size;
	stk = new Type[MAXN];
	top = -1;
}
template <typename Type>
Stack<Type>::~Stack() {
	delete stk;
}
template <typename Type>
int Stack<Type>::push(Type x) {
	if (top >= MAXN - 1)return -1;
	stk[++top] = x;
	return 0;
}
template <typename Type>
int Stack<Type>::pop(Type* x) {
	if (top == -1)return 0;
	*x = stk[top--];
	return 0;
}
template <typename Type>
Type Stack<Type>::getTop() {
	if (top == -1)return NULL;//null要大写
	return stk[top];
}
template <typename Type>
int Stack<Type>::isEmpty()const {
	return top == -1;
}
template <typename Type>
int Stack<Type>::isFull()const {
	return top == MAXN - 1;
}
template <typename Type>
int Stack<Type>::size() const {
	return top;
}
