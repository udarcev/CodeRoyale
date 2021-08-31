#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("omit-frame-pointer")
#pragma GCC optimize("unroll-loops")

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

high_resolution_clock::time_point now = high_resolution_clock::now();
#define TIME duration_cast<duration<double>>(high_resolution_clock::now() - now).count()

class Point;
class Site;
class Unit;
class Action;


void save();
void load();
void play(int);
double evaluate();
double evaluateExpand();
double evaluateSave();
double evaluateFarm();

bool needFarm = false;
bool needSave = false;
bool needBarracks = false;
bool panic = false;


inline int fastrand() {
    static unsigned int g_seed = 42;
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}

inline int rnd(int b) {
    return fastrand() % b;
}

inline int rnd(int a, int b) {  //both inclusive
    return a + rnd(b - a + 1);
}


constexpr double EPSILON = 0.00001;
constexpr int criticalTowerHealth = 400;
const uint8_t noStruct          = 0x01;
const uint8_t Goldmine          = 0x02;
const uint8_t Tower             = 0x04;
const uint8_t BarracksKnight    = 0x08;
const uint8_t BarracksArcher    = 0x10;
const uint8_t BarracksGiant     = 0x20;
const uint8_t option7           = 0x40;
const uint8_t option8           = 0x80;
const uint8_t anyBuilding = 255;
const uint8_t anyBarracks = BarracksGiant | BarracksArcher | BarracksKnight;
const uint8_t neutralOwner   = 0x01;
const uint8_t mineOwner      = 0x02;
const uint8_t enemyOwner     = 0x04;
const uint8_t anyOwner = 255;


enum Actions {NONE = -1, WAIT, MOVE, BUILD};
enum Buildings {MINE = 0, TOWER, KNIGHT_BARRACKS, ARCHER_BARRACKS, GIANT_BARRACKS};
enum UnitType {QUEEN = -1, KNIGHT, ARCHER, GIANT};
string unitsStr[] {"QUEEN", "KNIGHT", "ARCHER", "GIANT"};
string buildStr[] {"MINE", "TOWER", "BARRACKS-KNIGHT", "BARRACKS-ARCHER", "BARRACKS-GIANT"};


int gold;
int c_gold;
int touchedSite;

constexpr double WIDTH = 1920;
constexpr double HEIGHT = 1000;


vector<Site*> sites;
vector<Unit*> units;
Unit* myQueen = nullptr;
Unit* enQueen = nullptr;
Point* startCorner = nullptr;
Point* enemyCorner = nullptr;
Point* opStartCorner = nullptr;
Point* opEnemyCorner = nullptr;
Point* midStartCorner = nullptr;
Point* midEnemyCorner = nullptr;


constexpr int DEPTH = 8;
constexpr int POOL = 50;
constexpr int MUTATE_PROB = 35;
int sims = 0;
int all_sims = 0;

vector<Action*> possibleMoves;
Action* waitAction = nullptr;

vector<Site*> bestSitesPush;
vector<Site*> bestTowersRegen;
vector<Site*> bestMines;

const int dirCount = 16;                ///NEED TO EXPERIMENT WITH IT///
Point* g_points[dirCount];

int current_move = -1;

template<class T>
struct vec2 {
    typedef vec2<T> vec2d;
    T x, y;
    inline vec2<T>() : x(0), y(0) { }
    inline explicit vec2<T>(T v) : x(v), y(v) { }
    inline vec2<T>(T x, T y) : x(x), y(y) { }
    template<class TT> inline vec2<T>(vec2<TT> v) : x(T(v.x)), y(T(v.y)) {}
    inline vec2<T> rounded()const { return vec2<int>(round(x), round(y)); }
    inline void roundV() { x = round(x); y = round(y); }
    inline void clamp(int radius) {
        if (x < radius) x = radius; else if (x > WIDTH - radius) x = WIDTH - radius;
        if (y < radius) y = radius; else if (y > HEIGHT - radius) y = HEIGHT - radius;
    }
    inline T dot(const vec2d v)const { return x*v.x + y*v.y; }
    inline T cross(const vec2d v)const { return y*v.x - x*v.y; }
    inline T lengthSq()const { return dot(*this); }
    inline double length()const { return sqrt(lengthSq()); }
    inline vec2<T> normalized()const { return vec2<double>(x, y) / length(); }
    inline vec2<T> operator+(vec2d b)const { return vec2d(x + b.x, y + b.y); }
    inline vec2<T> operator-(vec2d b)const { return vec2d(x - b.x, y - b.y); }
    template<class T2> inline vec2<T> operator*(T2 n)const { return vec2d(x*n, y*n); }
    template<class T2> inline vec2<T> operator/(T2 n)const { return vec2d(x / n, y / n); }
    inline vec2<T> operator+=(vec2d b) { return *this = vec2d(x + b.x, y + b.y); }
    inline vec2<T> operator-=(vec2d b) { return *this = vec2d(x - b.x, y - b.y); }
    template<class T2> inline vec2<T> operator*=(T2 val) { return *this = vec2d(x * val, y * val); }
    template<class T2> inline vec2<T> operator/=(T2 val) { return *this = vec2d(x / val, y / val); }
    inline bool operator==(vec2d b) { return x == b.x && y == b.y; }
    inline bool operator!=(vec2d b) { return x != b.x || y != b.y; }
    inline void print() { cerr << x << " " << y << endl; }

    friend ostream &operator<<(ostream &os, const vec2 &vec21) {
        os << vec21.x << " " << vec21.y;
        return os;
    }
};
using vec2d = vec2<double>;
using vec2i = vec2<int>;


bool enemyKnightInRadius(Site* p);
Unit* findClosestUnit(Point* p, uint8_t owner, bool withoutQueen);
Site* findClosestSite(Point* p, uint8_t structureType, uint8_t owner);
int amountOfSites(uint8_t structureType, uint8_t owner = anyOwner);
Site* mineFromTower();
Site* towerToBuild = nullptr;
Site* barrackFromTower();
Site* barrackToBuild = nullptr;


class Point {
public:
    vec2d coord;

    Point() : coord(0, 0) {}
    Point(vec2d coord) : coord(coord) {}
    Point(double x, double y) : coord(x, y) {}

    inline virtual double dist(Point p) {
        return (coord - p.coord).length();
    }

    inline virtual double dist(Point *p) {
        return (coord - p->coord).length();
    }

    inline virtual double dist2(Point *p) {
        return (coord - p->coord).lengthSq();
    }

    inline virtual bool isInRange(Point* p, double range) {
        return (this->dist2(p)) <= (range * range);
    }

    virtual void moveAway(Point *p, double distance) {
        double d = this->dist2(p);
        if (d < EPSILON) {
            this->coord += vec2d(1., 0.) * (distance + 1);
            return;
        }
        this->coord -= (p->coord - coord).normalized() * (distance + 1);
    }

    bool operator==(const Point &rhs) const {
        return coord.x == rhs.coord.x && coord.y == rhs.coord.y;
    }
};


class Action {
public:
    int action;
    Point* coord;
    int siteId;
    int building;
    int health;

    Action(int action) : action(action) {}
    Action(int action, Point *coord) : action(action), coord(coord) {}
    Action(int action, int siteId, int building, int health = 780) : action(action), siteId(siteId), building(building), health(health) {}

    friend ostream &operator<<(ostream &os, const Action &action) {
        if (action.action == WAIT) {
            os << "WAIT";
        }
        else if (action.action == MOVE) {
            os << "MOVE " << (((Point*)myQueen)->coord + action.coord->coord).rounded();
        }
        else if (action.action == BUILD) {
            os << "BUILD " << action.siteId << " " << buildStr[action.building];
        }
        else {
            os << "WAIT";
        }
        return os;
    }

    bool operator==(const Action &rhs) const {
        return action == rhs.action &&
               *coord == *rhs.coord &&
               siteId == rhs.siteId &&
               building == rhs.building &&
               health == rhs.health;
    }
};


class Unit : public Point {
private:
    vec2d c_coord;
    int c_health;
    bool c_active;
    int c_respawn;
public:
    uint8_t owner;
    int unitType;
    int health;
    int radius;
    int moveSpeed;
    bool active = true;

    int respawn = -1;
    int barrackID = 0;

    Unit(double x, double y, uint8_t owner, int unitType, int health) : Point(x, y), owner(owner), unitType(unitType), health(health) {
        if (unitType == KNIGHT) {
            radius = 20;
            moveSpeed = 100;
        }
        else if (unitType == QUEEN) {
            radius = 30;
            moveSpeed = 60;
        }
        else if (unitType == ARCHER) {
            radius = 25;
            moveSpeed = 75;
        }
        else if (unitType == GIANT) {
            radius = 40;
            moveSpeed = 50;
        }
    }

    inline void save() {
        c_coord = coord;
        c_health = health;
        c_active = active;
        c_respawn = respawn;
    }

    inline void load() {
        coord = c_coord;
        health = c_health;
        active = c_active;
        respawn = c_respawn;
    }

    void prepareMoves(bool withMove = true, bool withTowers = false, bool withMines = false, bool withBarracks = false);

    void applyMove(bool bot);

    void applyMove(Action* action);

    void applyMove(Unit* u) {
        for (int i = 0; i < 5; i++) {
            int distance = (int)dist(u);
            if (distance <= u->radius + radius + 6)
                break;
            else {
                coord += (u->coord - coord).normalized() * min(moveSpeed / 5, distance - u->radius - radius);
                fixCollision();
            }
        }
        int distance = (int)dist(u);
        if (distance <= u->radius + radius + 6)
            u->health--;
    }

    void applyMoveAproximate(Unit* u) {
        int distance = (int)dist(u);
        if (distance > u->radius + radius + 5) {
            coord += (u->coord - coord).normalized() * min(moveSpeed, distance - u->radius - radius);
            fixCollision();
        }
    }

    void applyMove(Site* u);

    void breakSiteCheck();

    void play();

    inline void fixCollision();

    void end() {        // TODO округление координат
        if (respawn > -1) respawn--;
        if (active)
            coord.roundV();
        if (!active || unitType == QUEEN)
            return;
        health--;
        if (health <= 0)
            active = false;
    }

    friend ostream &operator<<(ostream &os, const Unit &unit) {
        os << "Type: " << unitsStr[unit.unitType + 1] << " coord: "  << unit.coord << " health: " << unit.health;
        return os;
    }
};


class Site : public Point {
private:
    int c_buildingType;
    uint8_t c_owner;
    uint8_t c_structureType;
    int c_goldRemaining;
    int c_param1;
    int c_param2;
public:
    int id;
    int radius;
    int goldRemaining = -1;
    int maxMineSize = -1;
    uint8_t structureType = noStruct;
    int buildingType = -1;
    uint8_t owner = neutralOwner;
    int param1 = -1;
    int param2 = -1;

    double calcWeight = 1;
    double towerWeight = 1;
    double mineWeight = 1;
    int timeToReach = 10;

    Site(double x, double y, int id, int radius) : Point(x, y), id(id), radius(radius) {}

    void update(int goldRemaining, int maxMineSize, uint8_t structureType, int buildingType, uint8_t owner, int param1, int param2) {
        if (this->goldRemaining != 0 && goldRemaining != -1)
            this->goldRemaining = goldRemaining;
        if (this->maxMineSize == -1)
            this->maxMineSize = maxMineSize;
        this->structureType = structureType;
        this->buildingType = buildingType;
        this->owner = owner;
        this->param1 = param1;
        this->param2 = param2;
        this->towerWeight = 1;
        this->mineWeight = 1;
    }

    void fixStruct() {
        switch (buildingType) {
            case NONE: structureType = noStruct; break;
            case MINE: structureType = Goldmine; break;
            case TOWER: structureType = Tower; break;
            case KNIGHT_BARRACKS: structureType = BarracksKnight; break;
            case ARCHER_BARRACKS: structureType = BarracksArcher; break;
            case GIANT_BARRACKS: structureType = BarracksGiant; break;
            default: structureType = noStruct;
        }
    }

    bool goodForTower() {
        if (owner & enemyOwner) {
            if (structureType & Tower)
                return false;
        }
        else if (owner & mineOwner) {
            if (buildingType == MINE && (goldRemaining > 10 || goldRemaining == -1) && !panic)
                return false;
            if ((structureType & anyBarracks)) {
                if (param1 > 0)
                    return false;
                if (!panic)
                    return false;
            }
        }
        return true;
    }

    bool goodForMine() {
        if (goldRemaining < 10 && goldRemaining != -1)
            return false;
        if (owner & enemyOwner) {
            if (structureType & Tower)
                return false;
        }
        else if (owner & mineOwner) {
            if ((structureType & anyBarracks))
                return false;
        }

        return true;
    }

    bool goodForBarrack() {
        if (owner & enemyOwner) {
            if (structureType & Tower)
                return false;
        }
        else if (owner & mineOwner) {
            if (buildingType == MINE  && amountOfSites(Goldmine, mineOwner) < 3)
                return false;
            if ((structureType & anyBarracks))
                return false;
        }

        return true;
    }

    void breakSite() {
        buildingType = -1;
        structureType = noStruct;
        owner = neutralOwner;
    }

    void buildSite(int b_type, uint8_t o) {
        if (buildingType != b_type || owner != o) {
            if (owner == o && (structureType & anyBarracks) && param1 > 0) {
                return;
            }
            buildingType = b_type;
            owner = o;
            fixStruct();
            if (buildingType == TOWER)
                param1 = 200;
            else if (buildingType == MINE)
                param1 = 1;
            else
                param1 = 0;
        }
        else {
            if (buildingType == TOWER) {
                param1 += 100;
                if (param1 > 800)
                    param1 = 800;
            }
            else if (buildingType == MINE) {
                param1++;
                if (param1 > maxMineSize)
                    param1 = maxMineSize;
            }
        }

    }

    inline void save() {
        c_buildingType = buildingType;
        c_owner = owner;
        c_structureType = structureType;
        c_goldRemaining = goldRemaining;
        c_param1 = param1;
        c_param2 = param2;
    }

    inline void load() {
        buildingType = c_buildingType;
        owner = c_owner;
        structureType = c_structureType;
        goldRemaining = c_goldRemaining;
        param1 = c_param1;
        param2 = c_param2;
    }

    void play() {
        if (buildingType == MINE && owner == mineOwner) {
            gold += param1;
            goldRemaining -= param1;
            if (goldRemaining <= 0) {
                gold += goldRemaining;
                goldRemaining = 0;
                breakSite();
            }
        }
        else if (buildingType == TOWER) {
            param2 = (int)sqrt((param1 * 1000 + M_PI * radius * radius) / M_PI);

            Unit* target = nullptr;
            if (owner == mineOwner) {
                target = findClosestUnit(this, enemyOwner, true);
                if (target == nullptr)
                    target = enQueen;
            }
            else {
                target = findClosestUnit(this, mineOwner, true);
                if (target == nullptr)
                    target = myQueen;
            }
            if (target != nullptr) {
                int distance = (int)dist(target);
                if (distance <= param2) {
                    int damage = 3 + (param2 - distance) / 200;
                    target->health -= damage;
                    if (target->unitType == QUEEN)
                        target->health += 2;
                }
                else if (target->unitType != QUEEN) {
                    target = (owner == mineOwner) ? enQueen : myQueen;
                    distance = (int)dist(target);
                    if (distance <= param2) {
                        int damage = 1 + (param2 - distance) / 200;
                        target->health -= damage;
                    }
                }
            }
        }
        else if ((structureType & anyBarracks) && param1 > 0) {
            param1--;
        }
    }

    void end() {
        if (buildingType == TOWER) {
            param1 -= 4;
            if (param1 <= 0) {
                breakSite();
            }
        }
    }

    friend ostream &operator<<(ostream &os, const Site &site) {
        os << " id: " << site.id << " goldRemaining: " << site.goldRemaining
           << " buildingType: " << buildStr[max(0, site.buildingType)] << " owner: " << ((site.owner & mineOwner) ? 0 : ((site.owner & enemyOwner) ? 1 : -1)) << " param1: " << site.param1
           << " param2: " << site.param2;
        return os;
    }
};

void Unit::breakSiteCheck() {
    Site* closest = findClosestSite(this, Goldmine, (owner == mineOwner ? enemyOwner : mineOwner));
    if (closest != nullptr && isInRange(closest, closest->radius + radius + 5))
        closest->breakSite();
}

void Unit::play()  {
    if (respawn == 0 && (sites[barrackID]->structureType & anyBarracks)) {
        active = true;
    }
    if (!active)
        return;
    if (unitType == KNIGHT) {        // TODO для своих тоже мб? сделать applyAprox для своих
        if (owner == enemyOwner) {
            Unit* target = myQueen;
            applyMove(target);
            breakSiteCheck();
        }
        else {
            //applyMoveAproximate(enQueen);
            applyMove(enQueen);
        }
    }
    else if (unitType == GIANT && owner == enemyOwner) {
        Site* target = findClosestSite(this, Tower, mineOwner);
        if (target != nullptr) {
            applyMove(target);
            breakSiteCheck();
        }
    }
    else if (unitType == ARCHER && owner == enemyOwner) {
        // По ходу игры решил полностью не симулировать их, так как это уменьшало количество симуляций, но предсказание этих ходов не давало преимущества
    }
    else if (unitType == QUEEN) {
        // По ходу игры решил полностью не симулировать их, так как это уменьшало количество симуляций, но предсказание этих ходов не давало преимущества
    }
}

inline void Unit::fixCollision() {
    Site* closest = findClosestSite(this, anyBuilding, anyOwner);
    if (closest == nullptr)
        return;

    double distance = dist(closest);
    if (distance < closest->radius + radius) {
        double delta = closest->radius + radius - distance;
        moveAway(closest, delta);
    }
}

void Unit::prepareMoves(bool withMove, bool withTowers, bool withMines, bool withBarracks) {
    possibleMoves.push_back(new Action(WAIT));
    if (withMove) {
        for (Point* p : g_points) {
            possibleMoves.push_back(new Action(MOVE, p));
        }
    }
    for (Site* s : sites) {
        if (dist(s) > 800)
            continue;

        if (withTowers) {       // TODO мб чекать что массив непустой и врубать панику?
            if (s->goodForTower())
                possibleMoves.push_back(new Action(BUILD, s->id, TOWER));
        }
    }
    if (withMines)
        possibleMoves.push_back(new Action(BUILD, towerToBuild->id, MINE));
    if (withBarracks)
        possibleMoves.push_back(new Action(BUILD, barrackToBuild->id, KNIGHT_BARRACKS));

    cerr << "possibleMoves: " << possibleMoves.size() << endl;
    // TODO сделать более умное добавление
}

void Unit::applyMove(bool bot) {
    Site* s = findClosestSite(this, anyBuilding, anyOwner);
    if (s != nullptr && ((s->buildingType == TOWER && (s->owner & owner)) || s->buildingType == -1 ) && isInRange(s, s->radius + radius + 5)) {
        sites[s->id]->buildSite(TOWER, owner);
    }
}

void Unit::applyMove(Action *action)  {
    if (action->action == WAIT)
        return;
    if (action->action == MOVE) {
        coord += action->coord->coord;

        Site* closest = findClosestSite(this, anyBarracks | Goldmine, enemyOwner);
        if (closest != nullptr && isInRange(closest, closest->radius + radius + 5))     // TODO по идее можно ломать прямо в фикс коллизии
            closest->breakSite();

        fixCollision();
        coord.clamp(radius);
    }
    else if (action->action == BUILD) {
        double distance = dist(sites[action->siteId]);
        if (distance < radius + sites[action->siteId]->radius + 5) {
            sites[action->siteId]->buildSite(action->building, owner);
        }
        else {
            coord += (sites[action->siteId]->coord - coord).normalized() * min(moveSpeed, (int)(dist(sites[action->siteId]) - sites[action->siteId]->radius - radius - 1));

            Site* closest = findClosestSite(this, anyBarracks | Goldmine, enemyOwner);
            if (closest != nullptr && isInRange(closest, closest->radius + radius + 5))     // TODO по идее можно ломать прямо в фикс коллизии
                closest->breakSite();

            fixCollision();
            //coord.clamp(radius);    // не может быть что он понадобится
        }
    }
}

void Unit::applyMove(Site *u) {
    for (int i = 0; i < 5; i++) {
        int distance = (int)dist(u);
        if (distance <= u->radius + radius + 6)
            break;
        else {
            coord += (u->coord - coord).normalized() * min(moveSpeed / 5, distance - u->radius - radius);
            fixCollision();
        }
    }
    int distance = (int)dist(u);
    if (distance <= u->radius + radius + 6) {
        u->param1 -= 80;
        if (u->param1 <= 0)
            u->breakSite();
    }
}

vector<Site*> findAllSites(uint8_t structureType, uint8_t owner = anyOwner) {
    vector<Site*> result;

    for (Site* s : sites) {
        if ((s->structureType & structureType) && (s->owner & owner)) {
            result.push_back(s);
        }
    }

    return result;
}

vector<Site*> findAllSitesInRadius(Point* p, double raduis, uint8_t structureType, uint8_t owner = anyOwner) {
    vector<Site*> result;

    for (Site* s : sites) {
        if ((s->structureType & structureType) && (s->owner & owner) && s->isInRange(p, raduis)) {
            result.push_back(s);
        }
    }

    return result;
}

Unit* findClosestUnit(Point* p, uint8_t owner, bool withoutQueen = false) {
    Unit* res = nullptr;
    double minDist = 10000000;
    for (Unit* u : units) {
        if (u->active && (u->owner & owner)) {
            if (withoutQueen && u->unitType == QUEEN)
                continue;
            double dist2 = p->dist2(u);
            if (dist2 < minDist) {
                minDist = dist2;
                res = u;
            }
        }
    }
    return res;
}

bool enemyKnightInRadius(Site* p) {
    for (Unit* u : units) {
        if (u->active && (u->owner & enemyOwner) && u->isInRange(p, p->radius + 20 + 5) ) {
            return true;
        }
    }
    return false;
}

int amountOfSites(uint8_t structureType, uint8_t owner) {
    return findAllSites(structureType, owner).size();
}

int countIncome(uint8_t own = mineOwner) {
    int res = 0;
    vector<Site*> mines = findAllSites(Goldmine, own);
    for (Site* s : mines) {
        res += s->param1;
    }
    return res;
}

int countUnits(int type, uint8_t owner) {
    int res = 0;
    for (Unit* u : units) {
        if (u->active && u->unitType == type && (u->owner & owner))
            res++;
    }
    return res;
}

vector<Unit*> getAllUnits(int type, uint8_t owner) {
    vector<Unit*> res;
    for (Unit* u : units) {
        if (u->active && u->unitType == type && (u->owner & owner))
            res.push_back(u);
    }
    return res;
}

vector<Unit*> getAllUnitsInRadius(int type, uint8_t owner, Point* p, double radius) {
    vector<Unit*> res;
    for (Unit* u : units) {
        if (u->active && u->unitType == type && (u->owner & owner) && u->isInRange(p, radius))
            res.push_back(u);
    }
    return res;
}

int countUnitsInRadius(int type, uint8_t owner, Point* target, double radius) {
    int res = 0;
    for (Unit* u : units) {
        if (u->active && u->unitType == type && (u->owner & owner) && u->isInRange(target, radius))
            res++;
    }
    return res;
}

int buildingBarracksPrediction(uint8_t owner) {
    int res = 0;
    for (Site* s : sites) {
        if ((s->structureType & BarracksKnight) && (s->owner & owner) && s->param1 > 0)
            res += 4;
    }
    return res;
}

Site* findClosestSite(Point* p, uint8_t structureType, uint8_t owner = anyOwner) {
    double minDist = 10000000;
    Site* result = nullptr;
    for (Site* s : sites) {
        if ( (s->structureType & structureType) && (s->owner & owner) ) {
            double dist = p->dist(s);
            if (dist < minDist) {
                minDist = dist;
                result = s;
            }
        }
    }
    return result;
}

Site* findBestMine() {      // TODO также искать среди башен, если она сзади и не нужна итд
    double minDist = 10000000;
    Site* result = nullptr;
    for (Site* s : bestMines) {
        if (s->goodForMine()) {
            double dist = myQueen->dist(s);
            if (dist < minDist) {
                minDist = dist;
                result = s;
            }
        }
    }
    return result;
}

Site* mineFromTower() {
    Site* res = nullptr;
    vector<Site*> mines = findAllSites(Goldmine, mineOwner);

    for (Site* s : mines) {     // TODO если их несколько вернуть ближайшую
        if (s->param1 < s->maxMineSize)
            return s;
    }

    vector<Site*> towers = findAllSites(Tower, mineOwner);
    if (towers.size() < 5)
        return res;

    vec2d coord;

    for (Site* s : towers) {
        coord += s->coord;
    }
    coord = coord / towers.size();

    Point midT(coord);
    vector<double> dist2c;
    vector<double> dist2t;
    vector<double> sumDist;
    Point targetCorner;
    Site* closestBar = findClosestSite(myQueen, anyBarracks, enemyOwner);
    if (closestBar == nullptr)
        targetCorner = *enemyCorner;
    else
        targetCorner = Point(closestBar->coord);

    for (Site* s : towers) {
        dist2c.push_back(s->dist(targetCorner));
        dist2t.push_back(s->dist(midT));
        s->calcWeight = (dist2c.back() + dist2t.back());
    }

    sort( towers.begin( ), towers.end( ), [ ]( const Site* lhs, const Site* rhs )
    {
        return lhs->calcWeight > rhs->calcWeight;
    });

    for (Site* s : towers) {
        if (s->goodForMine())
            return s;
    }

    return res;
}

Site* barrackFromTower() {
    Site* res = nullptr;
    vector<Site*> barr = findAllSites(anyBarracks, mineOwner);

    if (barr.size() > 1)
        return res;

    vector<Site*> towers = findAllSites(Tower, mineOwner);
    if (towers.size() < 5)
        return res;

    vec2d coord;

    for (Site* s : towers) {
        coord += s->coord;
    }
    coord = coord / towers.size();

    Point midT(coord);
    vector<double> dist2c;
    vector<double> dist2t;
    vector<double> sumDist;
    Point targetCorner = Point(enQueen->coord);

    for (Site* s : towers) {
        dist2c.push_back(s->dist(targetCorner));
        dist2t.push_back(s->dist(midT));
        s->calcWeight = dist2c.back();
    }

    sort( towers.begin( ), towers.end( ), [ ]( const Site* lhs, const Site* rhs )
    {
        return lhs->calcWeight < rhs->calcWeight;
    });

    for (Site* s : towers) {
        if (s->dist(midEnemyCorner) < 1200 && s->goodForBarrack())
            return s;
    }

    return res;
}

Site* findBestTower() {
    double minDist = 10000000;
    Site* result = nullptr;
    for (Site* s : sites) {
        if (s->buildingType == TOWER && (s->owner == enemyOwner || s->param1 > 300))
            continue;
        if (s->goodForTower()) {
            double dist = myQueen->dist(s);
            if (dist < minDist) {
                minDist = dist;
                result = s;
            }
        }
    }
    return result;
}

Site* findBestTowerRegenerate() {
    double minDist = 10000000, minDist2 = 10000000;
    Site* result = nullptr, *result2 = nullptr;
    for (Site* s : bestTowersRegen) {
        if (s->param1 <= criticalTowerHealth + 200) {        // TODO сделать иногда приоритет на почти мертвые тавера даже если они не ближайшие
            double dist = myQueen->dist(s);
            if (dist < minDist2) {
                minDist2 = dist;
                result2 = s;
            }
            if (s->param1 <= criticalTowerHealth && dist < minDist) {
                minDist = dist;
                result = s;
            }
        }
    }

    return result;
}


int timeToReachSite(Site* u) {
    int res = 10;
    for (int i = 0; i < 10; i++) {
        int distance = (int)myQueen->dist(u);
        if (distance <= u->radius + myQueen->radius + 5) {
            res = i;
            break;
        }
        else {
            myQueen->coord += (u->coord - myQueen->coord).normalized() * min(myQueen->moveSpeed, distance - u->radius - myQueen->radius);
            myQueen->fixCollision();
            myQueen->end();
        }
    }
    myQueen->load();
    return res;
}

vector<Site*> getBestSiteTowerPush() {
    vector<Site*> res;

    Point *startingPoint, *enemyPoint;
    startingPoint = myQueen;
    enemyPoint = enQueen;
    vector<Site*> allEBarracks = findAllSites(anyBarracks, enemyOwner);
    vector<Site*> allSitesB = findAllSitesInRadius(startingPoint, 800, anyBuilding, anyOwner);

    vector<Site*> allSites;
    for (Site* s : allSitesB) {
        if (s->buildingType != TOWER)
            allSites.push_back(s);
    }

    if (allSites.empty())
        return res;

    for (Site* s : allSites) {
        s->timeToReach = timeToReachSite(s);
        double tmpW = s->dist(enemyPoint);
        for (Site* other : allEBarracks)
            tmpW += s->dist(other);
        s->calcWeight = tmpW;
    }

    sort( allSites.begin( ), allSites.end( ), [ ]( const Site* lhs, const Site* rhs )
    {
        if (lhs->timeToReach == rhs->timeToReach)
            return lhs->calcWeight < rhs->calcWeight;
        return lhs->timeToReach < rhs->timeToReach;
    });

    vector<Site*> toRem;
    for (int i = allSites.size() - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (allSites[i]->calcWeight - allSites[j]->calcWeight > 200 || allSites[i]->timeToReach > 9) {
                toRem.push_back(allSites[i]);
                break;
            }
        }
    }

    for (Site* s : allSites) {
        bool push = true;
        for (Site* o : toRem) {
            if (s->id == o->id) {
                push = false;
                break;
            }
        }
        if (push)
            res.push_back(s);
    }

    double min = 100000;
    double max = 0;
    for (Site* s : res) {
        if (s->calcWeight > max)
            max = s->calcWeight;
        if (s->calcWeight < min)
            min = s->calcWeight;
    }

    for (Site* s : res) {
        s->calcWeight -= min;
        s->towerWeight = -s->calcWeight / (max - min) + 2;
    }

    return res;
}

vector<Site*> getBestTowerRegen() {
    vector<Site*> res;

    Point *startingPoint;
    startingPoint = myQueen;
    vector<Site*> allSites = findAllSitesInRadius(startingPoint, 800, Tower, mineOwner);

    if (allSites.empty())
        return res;

    for (Site* s : allSites) {
        s->timeToReach = timeToReachSite(s);
        double tmpW = 800 - s->param1;
        s->calcWeight = tmpW;
        if (tmpW >= 200)
            res.push_back(s);
    }

    return res;
}

vector<Site*> getBestMines() {
    vector<Site*> res;

    Point *startingPoint, *enemyPoint;
    startingPoint = myQueen;
    enemyPoint = enQueen;
    vector<Site*> allTowers = findAllSites(Tower, mineOwner);
    vector<Site*> allSitesB = findAllSitesInRadius(startingPoint, 800, anyBuilding, anyOwner);

    vector<Site*> allSites;
    for (Site* s : allSitesB) {
        if (s->buildingType == TOWER && s->owner == enemyOwner)
            continue;
        if (s->goldRemaining < 10 && s->goldRemaining != -1)
            continue;
        if (s->buildingType == MINE && s->owner == mineOwner && s->param1 >= s->maxMineSize)
            continue;

        allSites.push_back(s);
    }

    if (allSites.empty())
        return res;

    return res;
}

void predictMineSize() {
    for (int i = 0; i < sites.size(); i+=2) {
        if (sites[i]->maxMineSize != -1) {
            sites[i+1]->maxMineSize = sites[i]->maxMineSize;
        }
    }
    for (int i = 1; i < sites.size(); i+=2) {
        if (sites[i]->maxMineSize != -1) {
            sites[i-1]->maxMineSize = sites[i]->maxMineSize;
        }
    }
}

void save() {
    for (Site* s : sites) {
        s->save();
    }
    for (Unit* u : units) {
        u->save();
    }
    c_gold = gold;
}

void load() {
    for (Site* s : sites) {
        s->load();
    }
    for (Unit* u : units) {
        u->load();
    }
    gold = c_gold;
}

void play(int move = 0) {
    for (Unit* u : units) {
        u->play();
    }
    for (Site* s : sites) {
        s->play();
    }
    for (Unit* u : units) {
        u->end();
    }
    for (Site* s : sites) {
        s->end();
    }
    sims++;
}

bool shouldBuildTowers(int knights, int towers) {
    if (towers == 0)
        return true;
    if (towers <= 2) {
        if ((double)knights / 2 > towers)
            return true;
    }
    else if (towers < 6) {
        if ((double)knights / 2 > (towers + 1))
            return true;
    }

    return false;
}

// Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
double evaluate() {
    if (needFarm)
        return evaluateFarm();
    else if (needSave)
        return evaluateSave();
    else if (needBarracks)
        return evaluateExpand();

    return evaluateSave();
}

// Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
double evaluateSave() {
    double res = 0;

    for (Unit* u : units) {
        if (u->active && u->unitType == KNIGHT && (u->owner & enemyOwner)) {
            res -= u->health;
        }
    }

    return res;
}

// Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
double evaluateFarm() {
    double res = 0;

    vector<Site*> mines = findAllSites(Goldmine, mineOwner);
    for (Site* s : mines) {
        res += s->param1 * 100;
    }

    return res;
}

// Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
double evaluateExpand() {
    double res = 0;

    res += myQueen->health * 1000;

    return res;
}


/// Вот где то тут начинается всё для генетического алгоритма
class Solution {
public:
    static constexpr double minScore = -1000000;
    double score = Solution::minScore;
    Action* moves[DEPTH];

    Solution(bool with_rnd = false) {
        if (with_rnd) randomize();
    }

    void shift() {
        for (int i = 1; i < DEPTH; i++) {
            moves[i-1] = moves[i];
        }
        randomize(DEPTH-1, true);
        score = Solution::minScore;
    }

    inline void mutate() {
        randomize(rnd(DEPTH));
    }

    void mutate(Solution* child) {
        std::copy(begin(moves), end(moves), begin(child->moves));

        child->mutate();
        child->score = Solution::minScore;
    }

    void randomize(int idx, bool full = false) {
        moves[idx] = possibleMoves[rnd(possibleMoves.size())];
        score = Solution::minScore;
    }

    inline void randomize() {
        for (int i = 0; i < DEPTH; i++) randomize(i, true);
    }

    Solution* merge(Solution* solution) {
        Solution* child = new Solution();
        child->score = Solution::minScore;

        for (int i = 0; i < DEPTH; ++i) {
            if (rnd(2)) {
                child->moves[i] = solution->moves[i];
            } else {
                child->moves[i] = moves[i];
            }
        }

        return child;
    }

    void copy(Solution* solution) {
        for (int i = 0; i < DEPTH; ++i) {
            moves[i] = solution->moves[i];
        }

        this->score = solution->score;
    }
};

class Bot {
public:
    virtual void move(int turn) = 0;
};

class MonteBot : public Bot {
public:
    Solution sol;
    vector<Bot*> oppBots;
    int id;

    MonteBot(int id) {
        this->id = id;
    }

    void move(int turn) {
        move(&sol, turn);
    }

    void move(Solution *sol, int turn) {
        myQueen->applyMove(sol->moves[turn]);
        if (turn == 0)
            enQueen->applyMove(true);
    }

    // это решение на базе алгоритма монте карло
    void solve(double time, bool with_seed = false) {
        Solution best;
        if (with_seed) {
            best = sol;
            best.shift();
        } else {
            best.randomize();
        }
        get_score(&best);

        Solution child;
        while (TIME < time) {
            best.mutate(&child);
            for (int i = 1; i < DEPTH / 2; i++) {
                if (rnd(100) < MUTATE_PROB) child.mutate();
            }
            if (get_score(&child) > get_score(&best)) best = child;
        }
        sol.copy(&best);
    }

    // это решение на базе генетического алгоритма (слегка не по канону конечно, но так надо было)
    void solve_ga(double time, bool with_seed = false, bool keepON = false) {
        Solution best_old, *best;
        if (with_seed) {
            best_old.copy(&sol);
            if (!keepON)
                best_old.shift();
        } else {
            best_old.randomize();
        }
        get_score(&best_old);

        Solution *base;
        if (current_move) {
            base = new Solution();

            for (int j = 0; j < DEPTH; ++j) {
                base->moves[j] = best_old.moves[j];
            }
        }

        Solution** pool = new Solution*[POOL];
        Solution** newPool = new Solution*[POOL];
        Solution** temp;
        int counter = POOL;

        best = new Solution();
        Solution* sol = new Solution();
        sol->randomize();

        get_score(sol);
        pool[0] = sol;

        best->copy(sol);

        Solution* tempBest = sol;

        // First generation
        int startI = 1;

        if (current_move) {
            // Populate the POOL with some copy of the previous best one
            for (int i = startI; i < POOL / 4; ++i) {
                Solution* solution = new Solution();
                solution->copy(base);

                get_score(solution);

                if (solution->score > tempBest->score) {
                    tempBest = solution;
                }

                pool[i] = solution;
            }

            delete base;

            startI = POOL / 4;
        }

        for (int i = startI; i < POOL; ++i) {
            Solution* solution = new Solution();
            solution->randomize();

            get_score(solution);

            if (solution->score > tempBest->score) {
                tempBest = solution;
            }

            pool[i] = solution;
        }

        if (tempBest->score > best->score) {
            best->copy(tempBest);
        }
        tempBest = best;

        int generation = 1;
        int bestGeneration = 1;

        int poolFE = 0;
        while (TIME < time) {
            // New generation

            // Force the actual best with a mutation to be in the pool
            Solution* solution = new Solution();
            solution->copy(tempBest);
            solution->mutate();
            solution->mutate();
            get_score(solution);

            if (solution->score > tempBest->score) {
                tempBest = solution;
            }

            newPool[0] = solution;

            counter += 1;

            poolFE = 1;
            while (poolFE < POOL && (TIME < time)) {
                int aIndex = rnd(POOL);
                int bIndex;

                do {
                    bIndex = rnd(POOL);
                } while (bIndex == aIndex);

                int firstIndex = pool[aIndex]->score > pool[bIndex]->score ? aIndex : bIndex;

                do {
                    aIndex = rnd(POOL);
                } while (aIndex == firstIndex);

                do {
                    bIndex = rnd(POOL);
                } while (bIndex == aIndex || bIndex == firstIndex);

                int secondIndex = pool[aIndex]->score > pool[bIndex]->score ? aIndex : bIndex;

                Solution* child = pool[firstIndex]->merge(pool[secondIndex]);

                if (!rnd(3)) {
                    child->mutate();
                }

                get_score(child);

                if (child->score > tempBest->score) {
                    tempBest = child;
                }

                newPool[poolFE++] = child;

                counter += 1;
            }

            // Burn previous generation !!
            for (int i = 0; i < POOL; ++i) {
                delete pool[i];
            }

            temp = pool;
            pool = newPool;
            newPool = temp;

            if (tempBest->score > best->score) {
                best->copy(tempBest);
                bestGeneration = generation;
            }
            tempBest = best;

            generation += 1;
        }

        this->sol.copy(best);

        for (int i = 0; i < poolFE; ++i) {
            delete pool[i];
        }

        delete [] pool;
        delete [] newPool;
    }

    double get_score(Solution* sol) {
        if (sol->score == Solution::minScore) {
            sol->score = get_bot_score(sol);
        }

        return sol->score;
    }

    double get_bot_score(Solution* sol) {
        double score = 0;
        for (int i = 0; i < DEPTH; i++) {
            // утечки памяти и какие-то ненайденные ошибки мне сильно мешали, что пришлось даже сделать это
            if (    (sol->moves[i]->action == MOVE && sol->moves[i]->coord > (void*)0x7fff00000000) ||
                    (sol->moves[i]->action == BUILD && !(sol->moves[i]->siteId >= 0 && sol->moves[i]->siteId < 30) )    ) {
                cerr << "FUCK IT AGAIN" << endl;
                return Solution::minScore + 1;
            }
            move(sol, i);
            play(i);
            score += evaluate();
        }

        load();
        return score;
    }
};

// отдельная оптимизация, давшая мне огромный прирост в лидерборде за счет очень эффективного (нет) использования первой секунды времени игры
// это просчет супер оптимальных первых ходов для моей меты в игре
// из-за отсутствия времени это решение - дичайший копипаст, из-за которого тут насколько много строк кода
namespace pathFinder {
    constexpr int DEEP = 16;
    vector<Buildings> buildingType;
    int lastMove = 0;

    class Solution {
    public:
        static constexpr double minScore = -1000000;
        double score = Solution::minScore;
        Action* moves[DEEP];
        vector<int> passedID;

        Solution(bool with_rnd = false) {
            if (with_rnd) randomize();
        }

        void shift() {
            for (int i = 1; i < DEEP; i++) {
                moves[i-1] = moves[i];
            }
            randomize(DEEP-1, true);
            score = Solution::minScore;
            passedID.clear();
        }

        inline void mutate() {
            randomize(rnd(DEEP));
        }

        void mutate(Solution* child) {
            std::copy(begin(moves), end(moves), begin(child->moves));

            child->mutate();
            child->score = Solution::minScore;
            child->passedID.clear();
        }

        void randomize(int idx, bool full = false) {
            moves[idx] = possibleMoves[rnd(possibleMoves.size())];
            score = Solution::minScore;
            passedID.clear();
        }

        inline void randomize() {
            for (int i = 0; i < DEEP; i++) randomize(i, true);
        }

        Solution* merge(Solution* solution) {
            Solution* child = new Solution();
            child->score = Solution::minScore;
            passedID.clear();

            for (int i = 0; i < DEEP; ++i) {
                if (rnd(2)) {
                    child->moves[i] = solution->moves[i];
                } else {
                    child->moves[i] = moves[i];
                }
            }

            return child;
        }

        void copy(Solution* solution) {
            for (int i = 0; i < DEEP; ++i) {
                moves[i] = solution->moves[i];
            }
            passedID.clear();
            for (int i : solution->passedID)
                passedID.push_back(i);

            this->score = solution->score;
        }
    };

    class MonteBot {
    public:
        Solution sol;

        MonteBot() {}

        void move(int turn) {
            move(&sol, turn);
        }

        void move(Solution *sol, int turn) {
            myQueen->applyMove(sol->moves[turn]);
        }

        double eval(Solution* sol) {
            double res = 0;
            res += sol->passedID.size();
            return res;
        }

        // Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
        double evalFinal(Solution* sol, int turns) {
            double res = 0;
            res += sol->passedID.size();
            return res;
        }

        void solve_ga(double time, bool with_seed = false, bool keepON = false) {
            Solution best_old, *best;
            if (with_seed) {
                best_old.copy(&sol);
                if (!keepON)
                    best_old.shift();
            } else {
                best_old.randomize();
            }
            get_score(&best_old);

            Solution *base;
            if (keepON) {
                base = new Solution();
                base->copy(&best_old);
            }

            Solution** pool = new Solution*[POOL];
            Solution** newPool = new Solution*[POOL];
            Solution** temp;
            int counter = POOL;

            best = new Solution();
            Solution* sol = new Solution();
            sol->randomize();

            get_score(sol);
            pool[0] = sol;

            best->copy(sol);

            Solution* tempBest = sol;

            // First generation
            int startI = 1;

            if (keepON) {
                // Populate the POOL with some copy of the previous best one
                for (int i = startI; i < POOL / 4; ++i) {
                    Solution* solution = new Solution();
                    solution->copy(base);

                    get_score(solution);

                    if (solution->score > tempBest->score) {
                        tempBest = solution;
                    }

                    pool[i] = solution;
                }

                delete base;

                startI = POOL / 4;
            }

            for (int i = startI; i < POOL; ++i) {
                Solution* solution = new Solution();
                solution->randomize();

                get_score(solution);

                if (solution->score > tempBest->score) {
                    tempBest = solution;
                }

                pool[i] = solution;
            }

            if (tempBest->score > best->score) {
                best->copy(tempBest);
            }
            tempBest = best;

            int generation = 1;
            int bestGeneration = 1;

            int poolFE;
            while (TIME < time) {
                // New generation

                // Force the actual best with a mutation to be in the pool
                Solution* solution = new Solution();
                solution->copy(tempBest);
                solution->mutate();
                solution->mutate();
                get_score(solution);

                if (solution->score > tempBest->score) {
                    tempBest = solution;
                }

                newPool[0] = solution;

                counter += 1;

                poolFE = 1;
                while (poolFE < POOL && (TIME < time)) {
                    int aIndex = rnd(POOL);
                    int bIndex;

                    do {
                        bIndex = rnd(POOL);
                    } while (bIndex == aIndex);

                    int firstIndex = pool[aIndex]->score > pool[bIndex]->score ? aIndex : bIndex;

                    do {
                        aIndex = rnd(POOL);
                    } while (aIndex == firstIndex);

                    do {
                        bIndex = rnd(POOL);
                    } while (bIndex == aIndex || bIndex == firstIndex);

                    int secondIndex = pool[aIndex]->score > pool[bIndex]->score ? aIndex : bIndex;

                    Solution* child = pool[firstIndex]->merge(pool[secondIndex]);

                    if (!rnd(3)) {
                        child->mutate();
                    }

                    get_score(child);

                    if (child->score > tempBest->score) {
                        tempBest = child;
                    }

                    newPool[poolFE++] = child;

                    counter += 1;
                }

                // Burn previous generation !!
                for (int i = 0; i < POOL; ++i) {
                    delete pool[i];
                }

                temp = pool;
                pool = newPool;
                newPool = temp;

                if (tempBest->score > best->score) {
                    best->copy(tempBest);
                    bestGeneration = generation;
                }
                tempBest = best;

                generation += 1;
            }

            this->sol.copy(best);

            for (int i = 0; i < poolFE; ++i) {
                delete pool[i];
            }

            delete [] pool;
            delete [] newPool;
        }

        double get_score(Solution* sol) {
            if (sol->score == Solution::minScore) {
                sol->passedID.clear();
                sol->score = get_bot_score(sol);
            }

            return sol->score;
        }

        void temp(Solution* sol) {
            Site* s = findClosestSite(myQueen, anyBuilding, anyOwner);
            if ((int)myQueen->dist(s) <= s->radius + myQueen->radius + 3) {
                bool add = true;
                for (int i : sol->passedID) {
                    if (s->id == i) {
                        add = false;
                        break;
                    }
                }
                if (add) {
                    sol->passedID.push_back(s->id);
                }
            }
        }

        double get_bot_score(Solution* sol) {
            double score = 0;
            temp(sol);
            for (int i = 0; i < DEEP; i++) {
                move(sol, i);
                myQueen->end();
                temp(sol);

                score += eval(sol);
                if (sol->passedID.size() >= 5) {
                    lastMove = i;
                    score += evalFinal(sol, i);
                    myQueen->load();
                    return score;
                }
            }

            myQueen->load();
            return score;
        }
    };

    void assignRoles(Solution* sol) {
        vector<int> ids = sol->passedID;
        vector<double> distToMyCorner;
        vector<double> distToEnCorner;

        for (int i : ids) {
            distToMyCorner.push_back(sites[i]->dist(startCorner));
            distToEnCorner.push_back(sites[i]->dist(enemyCorner));
        }

        int barrackID = 2;
        double minD = 100000;
        for (int i = 0; i < distToEnCorner.size(); i++) {
            if (distToEnCorner[i] < minD) {
                minD = distToEnCorner[i];
                barrackID = i;
            }
        }

        buildingType.push_back(MINE);
        buildingType.push_back(MINE);
        buildingType.push_back(KNIGHT_BARRACKS);
        buildingType.push_back(TOWER);
        buildingType.push_back(TOWER);
    }
}

// Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
void detectSituation(int enemyKnights, int enemyKnightsInRadius, int towersInRad) {
    needSave = needFarm = needBarracks = panic = false;
    needFarm = true;
}

// Реальный код удален, чтобы при перезаливке моего решения меня не сместили с моего места, ну и правилами платформы это запрещено =)
void detectTemp() {
    needSave = needFarm = needBarracks = panic = false;
    needSave = true;
}


list<Action*> actToTest;
int firstTowerId, secondTowerID;
void getBaseList(pathFinder::MonteBot bot, int st, int en) {
    for (int i = st; i <= en; i++) {
        actToTest.push_back(bot.sol.moves[i]);
    }
    firstTowerId = bot.sol.passedID[3];
    secondTowerID = bot.sol.passedID[4];
    actToTest.push_back(new Action(BUILD, secondTowerID, TOWER));
    actToTest.push_back(new Action(BUILD, secondTowerID, TOWER));
    actToTest.push_back(new Action(BUILD, secondTowerID, TOWER));
    actToTest.push_back(new Action(BUILD, secondTowerID, TOWER));
    actToTest.push_back(new Action(BUILD, secondTowerID, TOWER));    // TODO попробовать WAIT
};

double specEval() {
    double res = 0;
    vector<Site*> allTowers = findAllSites(Tower, mineOwner);
    vector<Unit*> e_knights = getAllUnits(KNIGHT, enemyOwner);
    for (Unit* u : e_knights) {
        res -= u->health;
    }
    res -= e_knights.size() * 50;
    for (Site* s : allTowers) {
        res += s->param1 * 2;
    }
    res += myQueen->health * 1000;
    return res;
}


vector<double> allEvals;
void hpTEST() {

    double score = 0;
    for (Action* a : actToTest) {
        myQueen->applyMove(a);
        play(0);
    }
    allEvals.push_back(specEval());
    load();
}

int bestHealth = 250;
int bestAmount = 1;
void IMBASearch() {
    for (int i = 0; i < 6; i++) {
        hpTEST();
        actToTest.push_front(new Action(BUILD, firstTowerId, TOWER));
    }

    double biggest = -10000000;
    for (int i = 0; i < 6; i++) {
        if (allEvals[i] > biggest) {
            biggest = allEvals[i];
            bestAmount = i;
        }
    }
    bestHealth = 200 + bestAmount * 100 - 50;
    cerr << "IMBASearch " << bestAmount << " " << bestHealth << endl;
}


int sitesBuilt = 0;
int pathMove = 0;
int tmpID = 0;
bool dropped = true;
bool wasAttacked = true;
int atkCounter = 0;

int main()
{
    for (int i = 0; i < dirCount; i++) {
        double x = 60.0;
        double y = 0;
        double angle = (360.0 / dirCount * i) * M_PI / 180.0;

        double x1 = x * cos(angle);
        double y1 = x * sin(angle);
        g_points[i] = new Point(x1, y1);
        g_points[i]->coord.roundV();
    }

    waitAction = new Action(WAIT);

    MonteBot bot(0);
    pathFinder::MonteBot pathBot;

    int numSites;
    cin >> numSites; cin.ignore();
    sites.reserve(numSites);
    units.reserve(100);
    for (int i = 0; i < numSites; i++) {
        int siteId;
        int x;
        int y;
        int radius;
        cin >> siteId >> x >> y >> radius; cin.ignore();
        sites.push_back(new Site(x, y, siteId, radius));
    }

    // game loop
    while (numSites) {
        current_move++;
        cin >> gold >> touchedSite; cin.ignore();
        for (int i = 0; i < numSites; i++) {
            int siteId;
            int goldRemaining; // -1 if unknown
            int maxMineSize; // -1 if unknown
            int structureType; // -1 = No structure, 0 = Goldmine, 1 = Tower, 2 = Barracks
            int owner; // -1 = No structure, 0 = Friendly, 1 = Enemy
            int param1;
            int param2;
            cin >> siteId >> goldRemaining >> maxMineSize >> structureType >> owner >> param1 >> param2; cin.ignore();
            uint8_t flag = noStruct;
            uint8_t own = neutralOwner;
            int buildingType = -1;
            if (structureType == 0) {
                flag = Goldmine;
                buildingType = MINE;
            }
            else if (structureType == 1) {
                flag = Tower;
                buildingType = TOWER;
            }
            else if (structureType == 2) {
                if (param2 == 0) {
                    flag = BarracksKnight;
                    buildingType = KNIGHT_BARRACKS;
                    if (owner == 1 && param1 > 0) {
                        for (int j = 0; j < 4; j++) {
                            units.push_back(new Unit(sites[siteId]->coord.x, sites[siteId]->coord.y, enemyOwner, KNIGHT, 30));
                            units.back()->active = false;
                            units.back()->respawn = param1;
                            units.back()->barrackID = siteId;
                        }
                    }
                }
                else if (param2 == 1) {
                    flag = BarracksArcher;
                    buildingType = ARCHER_BARRACKS;
                }
                else if (param2 == 2) {
                    flag = BarracksGiant;
                    buildingType = GIANT_BARRACKS;      // TODO мб тоже спавнить?
                }
            }

            if (owner == 0)
                own = mineOwner;
            else if (owner == 1)
                own = enemyOwner;
            sites[siteId]->update(goldRemaining, maxMineSize, flag, buildingType, own, param1, param2);
        }
        int numUnits;
        cin >> numUnits; cin.ignore();
        for (int i = 0; i < numUnits; i++) {
            int x;
            int y;
            int owner;
            int unitType; // -1 = QUEEN, 0 = KNIGHT, 1 = ARCHER, 2 = GIANT
            int health;
            cin >> x >> y >> owner >> unitType >> health; cin.ignore();
            uint8_t own;
            if (owner == 0)
                own = mineOwner;
            else
                own = enemyOwner;
            units.push_back(new Unit(x, y, own, unitType, health));
            if (unitType == QUEEN) {
                if (owner == 0) {
                    myQueen = units.back();
                    if (startCorner == nullptr) {
                        if (myQueen->coord.length() < 1000) {
                            startCorner = new Point(0, 0);
                            opStartCorner = new Point(0, HEIGHT);
                            midStartCorner = new Point(0, HEIGHT / 2);
                            enemyCorner = new Point(WIDTH, HEIGHT);
                            opEnemyCorner = new Point(WIDTH, 0);
                            midEnemyCorner = new Point(WIDTH, HEIGHT / 2);
                        }
                        else {
                            startCorner = new Point(WIDTH, HEIGHT);
                            opStartCorner = new Point(WIDTH, 0);
                            midStartCorner = new Point(WIDTH, HEIGHT / 2);
                            enemyCorner = new Point(0, 0);
                            opEnemyCorner = new Point(0, HEIGHT);
                            midEnemyCorner = new Point(0, HEIGHT / 2);
                        }
                    }
                }
                else
                    enQueen = units.back();
            }
        }

        predictMineSize();
        save();
        now = high_resolution_clock::now();

        if (current_move == 0) {
            myQueen->prepareMoves();

            pathBot.solve_ga(0.25, false);
            pathBot.solve_ga(0.5, true, true);
            pathBot.solve_ga(0.85, true, true);

            for (int i : pathBot.sol.passedID)
                cerr << i << " ";
            cerr << endl << pathBot.sol.score << endl;

            pathBot.sol.score = Solution::minScore;
            pathBot.get_score(&pathBot.sol);

            pathFinder::assignRoles(&pathBot.sol);
            possibleMoves.clear();
        }

        bestSitesPush = getBestSiteTowerPush();
        bestTowersRegen = getBestTowerRegen();
        bestMines = getBestMines();

        double time_limit = 0.04;

        Action tmp = Action(NONE);
        if (current_move == 0 && !dropped) {
            tmpID = findClosestSite(myQueen, anyBuilding, anyOwner)->id;
        }
        if (!dropped) {
            if (sites[tmpID]->buildingType != KNIGHT_BARRACKS && current_move < 5)
                tmp = Action(BUILD, tmpID, KNIGHT_BARRACKS);
            else {
                dropped = true;
            }
        }

        if (sitesBuilt < 5) {
            Site* close = sites[pathBot.sol.passedID[sitesBuilt]];

            if ((int)myQueen->dist(close) <= myQueen->radius + close->radius + 4) {
                if (close->buildingType == -1) {
                    tmp = Action(BUILD, close->id, pathFinder::buildingType[sitesBuilt]);
                }
                else {
                    if (close->buildingType == MINE) {
                        if (close->param1 < close->maxMineSize)
                            tmp = Action(BUILD, close->id, pathFinder::buildingType[sitesBuilt]);
                        else
                            sitesBuilt++;
                    }
                    else if (close->buildingType == TOWER) {
                        if (sitesBuilt == 3 && actToTest.size() == 0) {
                            getBaseList(pathBot, pathMove, pathFinder::lastMove);
                            IMBASearch();
                        }
                        if ( (close->param1 < bestHealth && sitesBuilt == 3) || (close->param1 < 250 && sitesBuilt > 3))
                            tmp = Action(BUILD, close->id, pathFinder::buildingType[sitesBuilt]);
                        else
                            sitesBuilt++;
                    }
                    else
                        sitesBuilt++;
                }
            }
            if (tmp.action == NONE && sitesBuilt < 5 && pathMove < pathFinder::DEEP - 1)
                tmp = *pathBot.sol.moves[pathMove++];
        }

        int enemyKnights = buildingBarracksPrediction(enemyOwner) + countUnits(KNIGHT, enemyOwner);
        int enemyKnightsInRadius = countUnitsInRadius(KNIGHT, enemyOwner, myQueen, 600);

        int income = countIncome();

        vector<Site*> allMyBarracks = findAllSites(anyBarracks, mineOwner);
        int myBarracksCount = allMyBarracks.size();
        int enemyBarracksCount = amountOfSites(anyBarracks, enemyOwner);
        vector<Site*> allTowers = findAllSites(Tower, mineOwner);
        int towersInRad = 0;
        for (Site* s : allTowers) {
            if (myQueen->isInRange(s, s->param2)) towersInRad++;
        }
        int towerCount = allTowers.size();
        Site* bestTowerRegen = findBestTowerRegenerate();
        Site* bestTower = findBestTower();
        Site* bestMine = findBestMine();
        Site* closestFreeSite = findClosestSite(myQueen, noStruct, anyOwner);

        towerToBuild = mineFromTower();
        barrackToBuild = barrackFromTower();

        detectSituation(enemyKnights, enemyKnightsInRadius, towersInRad);

        if (wasAttacked) {
            needSave = needFarm = needBarracks = panic = false;    // TODO TEMP
            needSave = true;
            detectTemp();
        }

        if (needSave)           {myQueen->prepareMoves(true, true, false);                      cerr << "SAVING" << endl;}
        else if (needFarm)      {myQueen->prepareMoves(true, false, true);                      cerr << "FARMING" << endl;}
        else if (needBarracks)  {myQueen->prepareMoves(true, false, false, true);     cerr << "BARRACKS" << endl;}
        else                    {myQueen->prepareMoves(true, true, false);                      cerr << "SAVING" << endl;}

        bot.solve_ga(time_limit / 2, current_move > 0);
        int oldHP = myQueen->health;
        for (int i = 0; i < DEPTH; i++) {
            myQueen->applyMove(bot.sol.moves[i]);
            play(i);
        }
        double oldEval = evaluate();

        int deltaHP = oldHP - myQueen->health;
        int creeps = countUnitsInRadius(KNIGHT, enemyOwner, myQueen, myQueen->radius + 20 + 100);
        load();
        if (deltaHP > 4 || creeps > 2) {
            cerr << "PANIC MODE" << endl;
            panic = true;

            for (int i = 0; i < possibleMoves.size(); i++)
                delete (possibleMoves[i]);
            possibleMoves.clear();

            if (needSave)           {myQueen->prepareMoves(true, true, false);          cerr << "SAVING" << endl;}
            else if (needFarm)      {myQueen->prepareMoves(true, false, true);          cerr << "FARMING" << endl;}
            else if (needBarracks)  {myQueen->prepareMoves(true, false, false, true);   cerr << "BARRACKS" << endl;}
            else                    {myQueen->prepareMoves(true, true, false);          cerr << "SAVING" << endl;}
        }
        bot.solve_ga(time_limit, current_move > 0, true);

        if (tmp.action != NONE)
            cout << tmp << endl;
        else {
            Action* bestmove1 = bot.sol.moves[0];
            if (bestmove1->action == WAIT) {
                Site* s = findClosestSite(myQueen, anyBuilding, anyOwner);
                if (s != nullptr && myQueen->isInRange(s, myQueen->radius + s->radius + 5)) {
                    if (s->owner & mineOwner) {
                        if (s->buildingType == MINE)
                            cout << Action(BUILD, s->id, MINE) << endl;
                        else if (s->buildingType == TOWER)
                            cout << Action(BUILD, s->id, TOWER) << endl;
                        else
                            cout << *bestmove1 << endl;
                    }
                    else
                        cout << Action(BUILD, s->id, TOWER) << endl;
                }
                else
                    cout << *bestmove1 << endl;
            }
            else
                cout << *bestmove1 << endl;

            for (Action* a : bot.sol.moves) {
                cerr << *a << endl;
            }

            int tmph = myQueen->health;

            cerr << "myQueen->health " << myQueen->health;
            myQueen->applyMove(bestmove1);
            play();
            cerr << " " << myQueen->health << endl;

            if (tmph != myQueen->health)
                cerr << "PANIC!!!" << endl;

        }

        vector<Site*> brks = findAllSites(anyBarracks, mineOwner);
        for (Site* s : brks)
            s->calcWeight = s->dist(enQueen);
        sort( brks.begin( ), brks.end( ), [ ]( const Site* lhs, const Site* rhs )
        {
            return lhs->calcWeight < rhs->calcWeight;
        });

        if (gold >= 240)
            atkCounter = 0;

        if (brks.size() >= 1 && gold >= 80 && (brks[0]->calcWeight < 1000 || atkCounter < 3)) {
            cout << "TRAIN";

            if (brks.size() == 2) {
                if (gold >= 160)
                    for (Site* s : brks)
                        cout << " " << s->id;
            } else
                cout << " " << brks[0]->id;

            cout << endl;
            wasAttacked = true;
            atkCounter++;
        }
        else
            cout << "TRAIN" << endl;

        for (int i = 0; i < units.size(); i++) {
            delete(units[i]);
        }
        units.clear();
        for (int i = 0; i < possibleMoves.size(); i++)
            delete (possibleMoves[i]);
        possibleMoves.clear();

        cerr << "SIMS: " << sims << endl;
        if (current_move > 0) {
            all_sims += sims;
            cerr << "Average sims: " << all_sims / current_move << endl;
        }
        sims = 0;
        cerr << TIME << endl;
    }
}