#include <cstring>
#include <iostream>
using namespace std;

int GetLongestCommonSubString(const char *pStr1, const char *pStr2)
{
	/* 判断参数合法性 */
	if (pStr1 == NULL || pStr2 == NULL)
	{
		return -1;
	}

	int n = strlen(pStr1);
	int m = strlen(pStr2);
	int longestCommonSubString = 0;

	/* 申请辅助空间，并初始化为0 */
	int *LCS = new int[m];
	for (int i = 0; i < m; i++)
	{
		LCS[i] = 0;
	}

	/* 不断判断pStr[i] ?= pStr[j]，然后根据不同情况来更新LCS */
	for (int i = 0; i < n; i++)
	{
		for (int j = m - 1; j >= 0; j--)
		{
			if (pStr1[i] == pStr2[j])	/* 如果pStr1[i] == pStr2[j]，LCS[j] = LCS[j-1] + 1 */
			{
				if (j == 0)
				{
					LCS[j] = 1;
				}
				else
				{
					LCS[j] = LCS[j-1] + 1;
				}
			}
			else						/* 如果pStr1[i] != pStr2[j]，LCS[j] = 0 */
			{
				LCS[j] = 0;
			}

			/* 更新最长子串的长度 */
			if (LCS[j] > longestCommonSubString)
			{
				longestCommonSubString = LCS[j];
			}
		}
	}

	delete LCS;
	LCS = NULL;

	return longestCommonSubString;
}

void Test(const char *testName, const char *pStr1, const char *pStr2, int expectedLongestCommonSubString)
{
	cout << testName << " : ";
	if (GetLongestCommonSubString(pStr1, pStr2) == expectedLongestCommonSubString)
	{
		cout << "Passed." << endl;
	}
	else
	{
		cout << "Failed." << endl;
	}
}

int main()
{
	Test("Test1", "caba", "bab", 2);
	Test("Test2", "abcd", "efg", 0);
	Test("Test3", "abcde", "abcde", 5);
}


#include <cstdio>
#include <iostream>
using namespace std;

int max(int a, int b)
{
	return a > b ? a : b;
}

int GetLongestCommonSequence(const char *pStr1, const char *pStr2)
{
	/* 判断参数的合法性 */
	if (pStr1 == NULL || pStr2 == NULL)
	{
		return -1;
	}

	int m = strlen(pStr1);
	int n = strlen(pStr2);

	/* 申请二维空间LCS[m+1][n+1] */
	int **LCS = new int*[m+1];
	for (int i = 0; i < m + 1; i++)
	{
		LCS[i] = new int[n+1];
	}

	/* 分别对LCS[i][0], LCS[0][j]赋值为0 */
	for (int i = 0; i < m+1; i++)
	{
		LCS[i][0] = 0;
	}
	for (int j = 0; j < n+1; j++)
	{
		LCS[0][j] = 0;
	}

	/* 分别遍历两个字符串，并更新LCS[i][j] */
	for (int i = 1; i < m+1; i++)
	{
		for (int j = 1; j < n+1; j++)
		{
			if (pStr1[i-1] == pStr2[j-1])
			{
				LCS[i][j] = LCS[i-1][j-1] + 1;
			}
			else
			{
				LCS[i][j] = max(LCS[i-1][j], LCS[i][j-1]);
			}
		}
	}

	/* 获取最长公共子序列 */
	int longestCommonSequence = LCS[m][n];

	/* 删除动态空间 */
	for (int i = 0; i < m + 1; i++)
	{
		delete [] LCS[i];
		LCS[i] = NULL;
	}	
	delete []LCS;
	LCS = NULL;

	/* 返回最长公共子序列 */
	return longestCommonSequence;
}

void Test(const char *testName, const char *pStr1, const char *pStr2, int expectedLongestCommonSequence)
{
	cout << testName << " : ";
	if (GetLongestCommonSequence(pStr1, pStr2) == expectedLongestCommonSequence)
	{
		cout << "Passed." << endl;
	}
	else
	{
		cout << "Failed." << endl;
	}
}

int main()
{
	Test("Test1", "ABCBDAB", "BDCABA", 4);
	Test("Test2", "A", "A", 1);
	Test("Test3", "AB", "BC", 1);
}


#include <cstring>
#include <iostream>
using namespace std;

int min(int a, int b, int c)
{
	int min = a;
	if (min > b)
	{
		min = b;
	}
	if (min > c)
	{
		min = c;
	}
	return min;
}

int GetLeastestEditDistance(const char *pStr1, const char *pStr2)
{
	if (pStr1 == NULL || pStr2 == NULL)
	{
		return -1;
	}

	int m = strlen(pStr1);
	int n = strlen(pStr2);

	/* 申请动态空间LED[m+1][n+1] */
	int **LED = new int *[m+1];
	for (int i = 0; i < m+1; i++)
	{
		LED[i] = new int[n+1];
	}

	/* 赋值LED[i][0] = i， LED[0][j] = j */
	for (int i = 0; i < m+1; i++)
	{
		LED[i][0] = i;
	}
	for (int j = 0; j < n+1; j++)
	{
		LED[0][j] = j;
	}

	/* 计算LED[i][j] */
	for (int i = 1; i < m+1; i++)
	{
		for (int j = 1; j < n+1; j++)
		{
			if (pStr1[i-1] == pStr2[j-1])
			{
				LED[i][j] = min(LED[i-1][j-1], LED[i-1][j] + 1, LED[i][j-1] + 1);
			}
			else
			{
				LED[i][j] = min(LED[i-1][j-1]+1, LED[i-1][j] + 1, LED[i][j-1] + 1);
			}
		}
	}

	/* 获得最小的编辑距离 */
	int leastestEditDistance = LED[m][n];

	/* 释放动态空间 */
	for (int i = 0; i < m+1; i++)
	{
		delete [] LED[i];
		LED[i] = NULL;
	}
	delete []LED;
	LED = NULL;

	/* 返回最小编辑距离 */
	return leastestEditDistance;
}

void Test(const char *testName, const char *pStr1, const char *pStr2, int expectedLeastestEditDistance)
{
	cout << testName << " : ";
	if (GetLeastestEditDistance(pStr1, pStr2) == expectedLeastestEditDistance)
	{
		cout << "Passed." << endl;
	}
	else
	{
		cout << "Failed." << endl;
	}
}

int main()
{
	Test("Test1", "a", "b", 1);
	Test("Test2", "abdd", "aebdd", 1);
	Test("Test3", "travelling", "traveling", 1);
	Test("Test4", "abcd", "abcd", 0);
	Test("Test5", NULL, NULL, -1);
}