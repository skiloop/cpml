/* 
 * File:   Point.h
 * Author: skiloop
 *
 * Created on 2013年11月22日, 下午4:43
 */

#ifndef POINT_H
#define	POINT_H

class Point {
public:
    Point(unsigned i = 0, unsigned j = 0, unsigned k = 0);
    Point(const Point& orig);
    virtual ~Point();
    unsigned int x;
    unsigned int y;
    unsigned int z;

    /**
     * 
     * @param maxPoint
     * @return 
     */
    bool checkMax(const Point&maxPoint) const;

    /**
     * 
     * @param xMax
     * @param yMax
     * @param zMax
     * @return 
     */
    bool checkMax(unsigned xMax, unsigned yMax, unsigned zMax) const;

    void increaseX() {
        x++;
    }

    void increaseY() {
        y++;
    }

    void increaseZ() {
        z++;
    }

    void decreaseX() {
        x--;
    }

    void decreaseY() {
        y--;
    }

    void decreaseZ() {
        z--;
    }

    void operator=(Point const &obj);

    void setValue(unsigned i, unsigned j, unsigned k);

private:


};

#endif	/* POINT_H */

